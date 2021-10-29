import math, time, pynwb
import os.path

from hdmf.backends.hdf5.h5_utils import H5DataIO

from ReadIntanData import *
from ProcessData import *
from WriteNWB import *
from SetupResources import *
from ReadIntanHeader import *

def convert_to_nwb(intan_filename,
                   blocks_per_chunk=1000,
                   use_compression=True,
                   compression_level=4,
                   subject=None,
                   manual_start_time=None,
                   lowpass_description=None,
                   highpass_description=None,
                   merge=None
                  ):
    """ Convert the specified Intan file(s) to NWB format.
    
    Parameters
    ----------
    intan_filename : str
        Name of .rhd file to convert. If this is an 'info.rhd' file (not from the Traditional File Format), then
        other files in this directory with a .dat suffix will also be read as data sources
    blocks_per_chunk : int
        Number of data blocks that should be included in each chunk of data.
        Higher values require more RAM, but may be faster and more efficient
    use_compression : bool
        Whether data in written NWB file should be compressed. If so, 'compression_level' will determine the level of compression
    compression_level : int
        Int ranging from 0 to 9 indicating the level of 'gzip' compression.
        Higher values decrease written NWB file size, but may increase the amount of time required to convert.
    subject : pynwb.file.Subject or None
        If present, this subject object contains metadata about the subject from which this data was gathered.
        Not including this will result in an NWB file that is ineligible for the DANDI archive.
    manual_start_time : datetime.datetime or None
        If present, this contains the date and time that the recording session started.
        If not, an attempt will be made to parse the .rhd file name for a timestamp to use.
    lowpass_description : str or None
        If present, this describes the filter (type, order, cutoff frequency, etc.) used to generate lowpass data file.
        Only applies if lowpass data was saved ('one file per signal type' or 'one file per channel' file format).
    highpass_description : str or None
        If present, this describes the filter (type, order, cutoff frequency, etc.) used to generate highpass data file.
        Only applies if highpass data was saved ('one file per signal type' or 'one file per channel' file format).
    merge : bool or None
        If present, whether merging should be attempted with other Intan files in this directory.
    
    Returns
    -------
    None
    """
    # Start timing
    tic = time.time()
    
    # Open file
    filesize = os.path.getsize(intan_filename)
    
    # Read file header
    header = read_header(intan_filename)
    
    # If merging is desired, get list of other files that are mergeable (their headers are similar enough to allow it)
    if merge:
        mergeable_files = get_mergeable_files(header)
    
    # Output a summary of recorded data
    print_summary(header)
    
    # Calculate how many data blocks are present (assuming 'traditional' format - calculate for other formats later)
    bytes_per_block = get_bytes_per_data_block(header)
    fids = {}
    
    total_num_data_blocks, file_format = get_data_size(filesize, header, fids, bytes_per_block)
    
    # Calculate how many samples of each signal type are present
    total_num_amp_samples = header['num_samples_per_data_block'] * total_num_data_blocks
    
    if merge:
        total_num_amp_samples += merged_samples('amplifier', mergeable_files)
    
    # Get filename and attempt to get session start time from filename
    out_filename, session_start_time = parse_filename(intan_filename)
    
    # If manual start time was specified, overwrite the automatically generated start time
    if manual_start_time is not None:
        session_start_time = manual_start_time
    
    # Get notes from header
    if not (header['notes']['note1'] or header['notes']['note2'] or header['notes']['note3']):
        session_description = 'no description provided'
    else:
        session_description = header['notes']['note1'] + ', ' + header['notes']['note2'] + ', ' + header['notes']['note3']
    
    # Set up NWB file
    nwbfile = pynwb.NWBFile(session_description=session_description,
                            identifier=out_filename[:-4],
                            session_start_time=session_start_time,
                            subject=subject
                           )
    
    # If suitable, create an 'ecephys' Processing Module for low/highpass data
    if header['lowpass_present'] or header['highpass_present']:
        ecephys_module = nwbfile.create_processing_module(name='ecephys', description='software-filtered ecephys data')
    
    # Create 'device' object
    intan_device = create_intan_device(nwbfile, header)
    
    # Create 'electrode_table_region' object
    electrode_table_region = create_electrode_table_region(nwbfile, header, intan_device)
    
    # Initialize variables before conversion begins
    chunks_to_read = initialize_chunk_list(total_num_data_blocks, blocks_per_chunk)
    num_gaps = previous_timestamp = blocks_completed = 0
    previous_samples = [0] * header['num_amplifier_channels'] * 2
    
    rhd = header['filetype'] == 'rhd'
            
    # For each chunk in chunks_to_read, read the Intan data and write the NWB data
    for chunk in range(len(chunks_to_read)):
        
        # Number of data blocks in this chunk
        num_data_blocks = chunks_to_read[chunk]
        
        # Number of unique samples (per channel) in this chunk
        amp_samples_this_chunk = header['num_samples_per_data_block'] * num_data_blocks
        
        # Possible to do: could use something like 'original_time_series' in RHSutilities to recycle a single vector of timestamps
        # Pre-allocate memory for data
        data = preallocate_data(header, file_format, amp_samples_this_chunk)
        
        # Initialize indices used when looping through data blocks
        indices = initialize_indices(header['filetype'])
        
        # Read all blocks in this chunk
        for block in range(num_data_blocks):
            read_one_data_block(header, data, indices, fids, file_format)
            
        # Extract digital input/output channels to separate variables
        # Don't do this for One File Per Channel file format, because the data has already been separated by channel
        if file_format != 'per_channel':
            if header['num_board_dig_in_channels'] > 0:
                extract_digital_data(header, data['board_dig_in_raw'], data['board_dig_in_data'])
            if header['num_board_dig_out_channels'] > 0:
                extract_digital_data(header, data['board_dig_out_raw'], data['board_dig_out_data'])
            
        if not rhd:
            extract_stim_data(header, data)
        
        # Check for gaps in timestamps
        t_key= 't_amplifier' if rhd else 't'
        previous_timestamp, num_gaps = check_for_gaps(data[t_key], num_gaps, previous_timestamp, chunk)
        
        # Scale to SI units
        scale(header, data, file_format)
        
        # Process wideband data with a notch filter if appropriate
        wideband_filter_string, previous_samples = process_wideband(header, chunk, data, previous_samples)
        
        # Wrap data arrays
        wrapped_data = wrap_data_arrays(header=header,
                                        data=data,
                                        t_key=t_key,
                                        amp_samples_this_chunk=amp_samples_this_chunk,
                                        total_num_amp_samples=total_num_amp_samples,
                                        use_compression=use_compression,
                                        compression_level=compression_level)
        
        if chunk == 0:
          
            if header['num_amplifier_channels'] > 0:
                # Create ElectricalSeries for amplifier data
                amplifier_series = pynwb.ecephys.ElectricalSeries(name='ElectricalSeries',
                                           data=wrapped_data.data_amplifier,
                                           electrodes=electrode_table_region,
                                           filtering=wideband_filter_string,
                                           resolution=1.95e-7,
                                           timestamps=wrapped_data.t,
                                           comments='voltage data recorded from the amplifiers of an Intan Technologies chip',
                                           description='voltage data recorded from the amplifiers of an Intan Technologies chip')
                nwbfile.add_acquisition(amplifier_series)
                
                if not rhd:
                    if header['dc_amplifier_data_saved']:
                        # Create TimeSeries for dc amplifier data
                        dc_amplifier_series = pynwb.TimeSeries(name='TimeSeries_dc',
                                                         data=wrapped_data.data_dc_amplifier,
                                                         resolution=0.01923,
                                                         unit='volts',
                                                         timestamps=amplifier_series,
                                                         comments='DC electrical voltage data recorded from an Intan Technologies chip',
                                                         description='DC electrical voltage data recorded from an Intan Technologies chip')
                        nwbfile.add_acquisition(dc_amplifier_series)
                        
                    # Create TimeSeries for amp settle data
                    amp_settle_series = pynwb.TimeSeries(name='TimeSeries_amp_settle',
                                                   data=wrapped_data.data_amp_settle,
                                                   unit='digital event',
                                                   timestamps=amplifier_series,
                                                   comments='amplifier settle activity of an Intan Technologies chip',
                                                   description='amplifier settle activity of an Intan Technologies chip')
                    nwbfile.add_stimulus(amp_settle_series)
                    
                    # Create TimeSeries for charge recovery data
                    charge_recovery_series = pynwb.TimeSeries(name='TimeSeries_charge_recovery',
                                                        data=wrapped_data.data_charge_recovery,
                                                        unit='digital event',
                                                        timestamps=amplifier_series,
                                                        comments='charge recovery activity of an Intan Technologies chip',
                                                        description='charge recovery activity of an Intan Technologies chip')
                    nwbfile.add_stimulus(charge_recovery_series)
                    
                    # Create TimeSeries for compliance limit data
                    compliance_limit_series = pynwb.TimeSeries(name='TimeSeries_compliance_limit',
                                                         data=wrapped_data.data_compliance_limit,
                                                         unit='digital event',
                                                         timestamps=amplifier_series,
                                                         comments='compliance limit activity of an Intan Technologies chip',
                                                         description='compliance limit activity of an Intan Technologies chip')
                    nwbfile.add_stimulus(compliance_limit_series)
                    
                    # Create TimeSeries for stim data
                    stim_series = pynwb.TimeSeries(name='TimeSeries_stimulation',
                                             data=wrapped_data.data_stim,
                                             resolution=header['stim_step_size'],
                                             unit='amps',
                                             timestamps=amplifier_series,
                                             comments='current stimulation activity of an Intan Technologies chip',
                                             description='current stimulation activity of an Intan Technologies chip')
                    nwbfile.add_stimulus(stim_series)
                    
                
                if header['lowpass_present']:
                    # Create ElectricalSeries for lowpass data
                    lowpass_series = pynwb.ecephys.ElectricalSeries(name='ElectricalSeries_lowpass',
                                                                    data=wrapped_data.data_lowpass,
                                                                    electrodes=electrode_table_region,
                                                                    filtering=lowpass_description,
                                                                    resolution=1.95e-7,
                                                                    timestamps=wrapped_data.t_lowpass,
                                                                    comments='lowpass voltage data',
                                                                    description='lowpass voltage data')
                    nwbfile.processing['ecephys'].add(lowpass_series)
                    
                if header['highpass_present']:
                    # Create ElectricalSeries for highpass data
                    highpass_series=pynwb.ecephys.ElectricalSeries(name='ElectricalSeries_highpass',
                                                                   data=wrapped_data.data_highpass,
                                                                   electrodes=electrode_table_region,
                                                                   filtering=highpass_description,
                                                                   resolution=1.95e-7,
                                                                   timestamps=wrapped_data.t,
                                                                   comments='highpass voltage data',
                                                                   description='highpass voltage data')
                    nwbfile.processing['ecephys'].add(highpass_series)
            
            if rhd:
                if header['num_aux_input_channels'] > 0:
                    # Create TimeSeries for auxiliary input data
                    aux_input_series = pynwb.TimeSeries(name='TimeSeries_aux_input',
                                                  data=wrapped_data.data_aux_in,
                                                  resolution=37.4e-6,
                                                  unit='volts',
                                                  timestamps=wrapped_data.t_aux_input,
                                                  comments='voltage data recorded from the auxiliary inputx of an Intan Technologies chip',
                                                  description='voltage data recorded from the auxiliary input of an Intan Technologies chip')
                    nwbfile.add_acquisition(aux_input_series)      

                if header['num_supply_voltage_channels'] > 0:
                    # Create TimeSeries for supply voltage data
                    supply_voltage_series = pynwb.TimeSeries(name='TimeSeries_supply_voltage',
                                                       data=wrapped_data.data_supply_voltage,
                                                       resolution=74.8e-6,
                                                       unit='volts',
                                                       timestamps=wrapped_data.t_supply_voltage,
                                                       comments='supply voltage data recorded from an Intan Technologies chip',
                                                       description='supply voltage data recorded from an Intan Technologies chip')
                    nwbfile.add_acquisition(supply_voltage_series)
                
            if header['num_board_adc_channels'] > 0:
                
                # Determine resolution for board adc TimeSeries
                if rhd:
                    if header['board_mode'] == 1:
                        resolution = 152.59e-6
                    elif header['board_mode'] == 13:
                        resolution = 312.5e-6
                    else:
                        resolution = 50.354e-6
                else:
                    resolution=312.5e-6
                    
                # If the amplifier ElectricalSeries has already been created, recycle that for its timestamps
                if header['num_amplifier_channels'] > 0:
                    board_adc_timestamps = amplifier_series
                    
                # Otherwise, set up board adc timestamps
                else:
                    board_adc_timestamps = wrapped_data.t
                
                # Create TimeSeries for board adc data
                board_adc_series = pynwb.TimeSeries(name='TimeSeries_analog_input',
                                              data=wrapped_data.data_board_adc,
                                              resolution=resolution,
                                              unit='volts',
                                              timestamps=board_adc_timestamps,
                                              comments='analog input data recorded from an Intan Technologies system',
                                              description='analog input data recorded from an Intan Technologies system')
                nwbfile.add_acquisition(board_adc_series)
                
            if not rhd:
                if header['num_board_dac_channels'] > 0:
                    
                    if header['num_amplifier_channels'] > 0:
                        board_dac_timestamps = amplifier_series
                        
                    else:
                        board_dac_timestamps = wrapped_data.t
                    
                    board_dac_series = pynwb.TimeSeries(name='TimeSeries_analog_output',
                                                  data=wrapped_data.data_board_dac,
                                                  resolution=312.5e-6,
                                                  unit='volts',
                                                  timestamps=board_dac_timestamps,
                                                  comments='analog output data recorded from an Intan Technologies system',
                                                  description='analog output data recorded from an Intan Technologies system')
                    nwbfile.add_acquisition(board_dac_series)
                
            if header['num_board_dig_in_channels'] > 0:
                
                # If the amplifier ElectricalSeries has already been created, recycle that for its timestamps
                if header['num_amplifier_channels'] > 0:
                    board_dig_in_timestamps = amplifier_series
                    
                # Otherwise,
                else:
                    # If the board adc TimeSeries has already been created, recycle that for its timestamps
                    if header['num_board_adc_channels'] > 0:
                        board_dig_in_timestamps = board_adc_series
                        
                    # Otherwise, set up board dig in timestamps
                    else:
                        board_dig_in_timestamps = wrapped_data.t
                
                # Create TimeSeries for digital input data
                board_dig_in_series = pynwb.TimeSeries(name='TimeSeries_digital_input',
                                                 data=wrapped_data.data_board_dig_in,
                                                 unit='digital event',
                                                 timestamps=board_dig_in_timestamps,
                                                 comments='digital input data recorded from an Intan Technologies system',
                                                 description='digital input data recorded from an Intan Technologies system')
                nwbfile.add_acquisition(board_dig_in_series)
                
            if header['num_board_dig_out_channels'] > 0:
                
                # If the amplifier ElectricalSeries has already been created, recycle that for its timestamps
                if header['num_amplifier_channels'] > 0:
                    board_dig_out_timestamps = amplifier_series
                    
                # Otherwise,
                else:
                    # If the board adc TimeSeries has already been created, recycle that for its timestamps
                    if header['num_board_adc_channels'] > 0:
                        board_dig_out_timestamps = board_adc_series
                        
                    # Otherwise,
                    else:
                        # If the board dig in TimeSeries has already been created, recycle that for its timestamps
                        if header['num_board_dig_in_channels'] > 0:
                            board_dig_out_timestamps = board_dig_in_series
                        
                        # Otherwise, set up board dig out timestamps
                        else:
                            board_dig_out_timestamps = wrapped_data.t
                
                # Create TimeSeries for digital output data
                board_dig_out_series = pynwb.TimeSeries(name='TimeSeries_digital_output',
                                                  data=wrapped_data.data_board_dig_out,
                                                  unit='digital event',
                                                  timestamps=board_dig_out_timestamps,
                                                  comments='digital output data recorded from an Intan Technologies system',
                                                  description='digital output data recorded from an Intan Technologies system')
                nwbfile.add_acquisition(board_dig_out_series)
                
            if rhd:
                if header['num_temp_sensor_channels'] > 0:

                    # If the supply voltage TimeSeries has already been created, recycle that for its timestamps
                    if header['num_supply_voltage_channels'] > 0:
                        temp_sensor_timestamps = supply_voltage_series

                    # Otherwise, use temp sensor timestamps
                    else:
                        temp_sensor_timestamps = wrapped_data.t_supply_voltage

                    # Create TimeSeries for temp sensor data
                    temp_sensor_series = pynwb.TimeSeries(name='TimeSeries_temperature_sensor',
                                                    data=wrapped_data.data_temp,
                                                    unit='deg C',
                                                    timestamps=temp_sensor_timestamps,
                                                    comments='temperature sensor data recorded from an Intan Technologies chip',
                                                    description='temperature sensor data recorded from an Intan Technologies chip')
                    nwbfile.add_acquisition(temp_sensor_series)
            
            # Write the data to file
            io = pynwb.NWBHDF5IO(out_filename, 'w')
            io.write(nwbfile)
            io.close()
            
        else:
            with pynwb.NWBHDF5IO(out_filename, mode='a') as io:
                append_nwbfile = io.read()
                
                # Append amplifier data
                if header['num_amplifier_channels'] > 0:
                    append_to_dataset(append_nwbfile.acquisition['ElectricalSeries'].timestamps, wrapped_data.t)
                    append_to_dataset(append_nwbfile.acquisition['ElectricalSeries'].data, wrapped_data.data_amplifier)
                    
                    if header['lowpass_present']:
                        append_to_dataset(append_nwbfile.processing['ecephys']['ElectricalSeries_lowpass'].timestamps, wrapped_data.t_lowpass)
                        append_to_dataset(append_nwbfile.processing['ecephys']['ElectricalSeries_lowpass'].data, wrapped_data.data_lowpass)
                        
                    if header['highpass_present']:
                        append_to_dataset(append_nwbfile.processing['ecephys']['ElectricalSeries_highpass'].timestamps, wrapped_data.t)
                        append_to_dataset(append_nwbfile.processing['ecephys']['ElectricalSeries_highpass'].data, wrapped_data.data_highpass)
                        
                    if not rhd:
                        if header['dc_amplifier_data_saved']:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_dc'].data, wrapped_data.data_dc_amplifier)
                            
                        append_to_dataset(append_nwbfile.stimulus['TimeSeries_amp_settle'].data, wrapped_data.data_amp_settle)
                        append_to_dataset(append_nwbfile.stimulus['TimeSeries_charge_recovery'].data, wrapped_data.data_charge_recovery)
                        append_to_dataset(append_nwbfile.stimulus['TimeSeries_compliance_limit'].data, wrapped_data.data_compliance_limit)
                        append_to_dataset(append_nwbfile.stimulus['TimeSeries_stimulation'].data, wrapped_data.data_stim)
                    
                if rhd:
                    # Append aux input data
                    if header['num_aux_input_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_aux_input'].timestamps, wrapped_data.t_aux_input)
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_aux_input'].data, wrapped_data.data_aux_in)

                    # Append supply voltage data
                    if header['num_supply_voltage_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_supply_voltage'].timestamps, wrapped_data.t_supply_voltage)
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_supply_voltage'].data, wrapped_data.data_supply_voltage)
                        
                    # Append temp sensor data
                    if header['num_temp_sensor_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_temperature_sensor'].data, wrapped_data.data_temp)
                        # If the timestamps vector hasn't been already appended via supply voltage data, append it here
                        if header['num_supply_voltage_channels'] == 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_temperature_sensor'].timestamps, wrapped_data.t_supply_voltage)
                            
                else:
                    # Append board dac data
                    if header['num_board_dac_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_output'].data, wrapped_data.data_board_dac)
                        # If the timestamps vector hasn't already been appended via amplifier data, append it here
                        if header['num_amplifier_channels'] == 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_output'].timestamps, wrapped_data.t)
                    
                # Append board adc data
                if header['num_board_adc_channels'] > 0:
                    append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_input'].data, wrapped_data.data_board_adc)
                    # If the timestamps vector hasn't been already appended via amplifier data, append it here
                    if header['num_amplifier_channels'] == 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_input'].timestamps, wrapped_data.t)
                    
                # Append board dig in data
                if header['num_board_dig_in_channels'] > 0:
                    append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_input'].data, wrapped_data.data_board_dig_in)
                    # If the timestamps vector hasn't been already appended via amplifier data or adc data, append it here
                    if header['num_amplifier_channels'] == 0:
                        if header['num_board_adc_channels'] == 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_input'].timestamps, wrapped_data.t)
                    
                # Append board dig out data
                if header['num_board_dig_out_channels'] > 0:
                    append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_output'].data, wrapped_data.data_board_dig_out)
                    # If the timestamps vector hasn't been already appended via amplifier data, adc data, or dig in data, append it here
                    if header['num_amplifier_channels'] == 0:
                        if header['num_board_adc_channels'] == 0:
                            if header['num_board_dig_in_channels'] == 0:
                                append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_output'].timestamps, wrapped_data.t)
                
                io.write(append_nwbfile)                 
        
        blocks_completed = blocks_completed + num_data_blocks
        percent_done = (blocks_completed / total_num_data_blocks) * 100
        print('Completed chunk {} ... {:0.2f}% done...'.format(chunk, percent_done))
        
    # Report whether gaps in timestamp data were found.
    if num_gaps == 0:
        print('No missing timestamps in data.')
    else:
        print('Warning: {} gaps in timestamp data found. Time scale will not be uniform!'.format(num_gaps))
            
    # Make sure we have read exactly the right amount of data
    bytes_remaining = filesize - header['fid'].tell()
    if bytes_remaining != 0:
        raise Exception('Error: End of file not reached.')
    
    # Close Intan file(s)
    for fid in fids.values():
        if type(fid) != bool: # aux_in_amplifier is a boolean value in the fids dictionary, so don't treat it as a fid
            fid.close()
        
    if merge:
        # Possible to do: this section is basically an abridged repeat of original file reading,
        # so it may be worth creating some high-level functions to allow for code reuse
        for header in mergeable_files:
            filesize = os.path.getsize(header['filename'])
            
            # Output a message announcing merging of this file
            print('\nMerging file: {}'.format(header['filename']))
            
            # Calculate how many data blocks are present (assuming 'traditional' format)
            bytes_per_block = get_bytes_per_data_block(header)
            fids = {}
            total_num_data_blocks, file_format = get_data_size(filesize, header, fids, bytes_per_block, False)
            
            # Only traditional files are mergeable
            if file_format != 'traditional':
                continue
                
            # Initialize variables before conversion begins
            chunks_to_read = initialize_chunk_list(total_num_data_blocks, blocks_per_chunk)
            num_gaps = previous_timestamp = blocks_completed = 0
            previous_samples = [0] * header['num_amplifier_channels'] * 2
            
            # For each chunk in chunks_to_read, read the Intan data and write the NWB data
            for chunk in range(len(chunks_to_read)):
                
                # Number of data blocks in this chunk
                num_data_blocks = chunks_to_read[chunk]
                
                # Number of unique samples (per channel) in this chunk
                amp_samples_this_chunk = header['num_samples_per_data_block'] * num_data_blocks
                
                # Pre-allocate memory for data
                data = preallocate_data(header, file_format, amp_samples_this_chunk)
                
                # Initialize indices used when looping through data blocks
                indices = initialize_indices(header['filetype'])
                
                # Read all blocks in this chunk
                for block in range(num_data_blocks):
                    read_one_data_block(header, data, indices, fids, file_format)
                    
                # Extract digital input/output channels to separate variables
                if header['num_board_dig_in_channels'] > 0:
                    extract_digital_data(header, data['board_dig_in_raw'], data['board_dig_in_data'])
                if header['num_board_dig_out_channels'] > 0:
                    extract_digital_data(header, data['board_dig_out_raw'], data['board_dig_out_data'])
                
                if not rhd:
                    extract_stim_data(header, data)
                
                # Check for gaps in timestamps
                t_key = 't_amplifier' if rhd else 't'
                previous_timestamp, num_gaps = check_for_gaps(data[t_key], num_gaps, previous_timestamp, chunk)
                
                # Scale to SI units
                scale(header, data, file_format)
                
                # Process wideband data with a notch filter if appropriate
                wideband_filter_string, previous_samples = process_wideband(header, chunk, data, previous_samples)
                
                # Wrap data arrays
                wrapped_data = wrap_data_arrays(header=header,
                                                data=data,
                                                t_key=t_key,
                                                amp_samples_this_chunk=amp_samples_this_chunk,
                                                total_num_amp_samples=total_num_amp_samples,
                                                use_compression=use_compression,
                                                compression_level=compression_level)
                    
                with pynwb.NWBHDF5IO(out_filename, mode='a') as io:
                    append_nwbfile = io.read()
                    
                    # Append amplifier data
                    if header['num_amplifier_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['ElectricalSeries'].timestamps, wrapped_data.t)
                        append_to_dataset(append_nwbfile.acquisition['ElectricalSeries'].data, wrapped_data.data_amplifier)
                        
                        # Don't bother with low or highpass, because merging necessitates traditional file format
                        if not rhd:
                            if header['dc_amplifier_data_saved']:
                                append_to_dataset(append_nwbfile.acquisition['TimeSeries_dc'].data, wrapped_data.data_dc_amplifier)
                                
                            append_to_dataset(append_nwbfile.stimulus['TimeSeries_amp_settle'].data, wrapped_data.data_amp_settle)
                            append_to_dataset(append_nwbfile.stimulus['TimeSeries_charge_recovery'].data, wrapped_data.data_charge_recovery)
                            append_to_dataset(append_nwbfile.stimulus['TimeSeries_compliance_limit'].data, wrapped_data.data_compliance_limit)
                            append_to_dataset(append_nwbfile.stimulus['TimeSeries_stimulation'].data, wrapped_data.data_stim)
                        
                    if rhd:
                        # Append aux input data
                        if header['num_aux_input_channels'] > 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_aux_input'].timestamps, wrapped_data.t_aux_input)
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_aux_input'].data, wrapped_data.data_aux_in)

                        # Append supply voltage data
                        if header['num_supply_voltage_channels'] > 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_supply_voltage'].timestamps, wrapped_data.t_supply_voltage)
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_supply_voltage'].data, wrapped_data.data_supply_voltage)
                            
                        # Append temp sensor data
                        if header['num_temp_sensor_channels'] > 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_temperature_sensor'].data, wrapped_data.data_temp)
                            # If the timestamps vector hasn't been already appended via supply voltage data, append it here
                            if header['num_supply_voltage_channels'] == 0:
                                append_to_dataset(append_nwbfile.acquisition['TimeSeries_temperature_sensor'].timestamps, wrapped_data.t_supply_voltage)
                                
                    else:
                        # Append board dac data
                        if header['num_board_dac_channels'] > 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_output'].data, wrapped_data.data_board_dac)
                            # If the timestamps vector hasn't already been appended via amplifier data, append it here
                            if header['num_amplifier_channels'] == 0:
                                append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_output'].timestamps, wrapped_data.t)
                        
                    # Append board adc data
                    if header['num_board_adc_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_input'].data, wrapped_data.data_board_adc)
                        # If the timestamps vector hasn't been alread appended via amplifier data, append it here
                        if header['num_amplifier_channels'] == 0:
                            append_to_dataset(append_nwbfile.acquisition['TimeSeries_analog_input'].timestamps, wrapped_data.t)
                            
                    # Append board dig in data
                    if header['num_board_dig_in_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_input'].data, wrapped_data.data_board_dig_in)
                        # If the timestamps vector hasn't been already appended via amplifier data or adc data, append it here
                        if header['num_amplifier_channels'] == 0:
                            if header['num_board_adc_channels'] == 0:
                                append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_input'].timestamps, wrapped_data.t)
                                
                    # Append board dig out data
                    if header['num_board_dig_out_channels'] > 0:
                        append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_output'].data, wrapped_data.data_board_dig_out)
                        # If the timestamps vector hasn't been already appended via amplifier data, adc data, or dig in data, append it here
                        if header['num_amplifier_channels'] == 0:
                            if header['num_board_adc_channels'] == 0:
                                if header['num_board_dig_in_channels'] == 0:
                                    append_to_dataset(append_nwbfile.acquisition['TimeSeries_digital_output'].timestamps, wrapped_data.t)
                                    
                        
                    io.write(append_nwbfile)
                
                blocks_completed = blocks_completed + num_data_blocks
                percent_done = (blocks_completed / total_num_data_blocks) * 100
                print('Completed chunk {} ... {:0.2f}% done...'.format(chunk, percent_done))
            
    
    print('Done! Elapsed time: {:0.2f} seconds'.format(time.time() - tic))