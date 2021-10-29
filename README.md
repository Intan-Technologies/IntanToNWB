# IntanToNWB
Convert Intan (rhd or rhs) files to NWB file format

## Installation
To run this example, either download the 8 files:
1. IntanToNWB.ipynb
2. ConverterUI.py
3. ConvertIntanToNWB.py
4. ProcessData.py
5. ReadIntanData.py
6. ReadIntanHeader.py
7. SetupResources.py
8. WriteNWB.py

to the same directory and open them with Jupyter Notebook on your own machine (Jupyter must already be installed), or click the Binder link below to run it in a remote Docker container.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/adrian-foy/IntanToNWB/HEAD)

This notebook depends on PyNWB, a Python NWB package that provides the HDF5 interface used for writing to NWB format. This notebook was written using version 2.0.0 of PyNWB.

## Usage
To read data from Intan files (.rhd or .rhs, or for other file formats, .dat), bring those files to the same directory this notebook is in, and change the "File to convert" to the file you wish to convert. After a period of time (progress can be monitored by reading the output print statements), an output file with the extension .nwb will be written in this directory.

If you wish to merge multiple continuous files into a single .NWB file, give the name of the earliest file in that recording session, and check "Merge Multiple Continuous Files" in the Advanced Settings section, and any files that continue from that recording session will automatically be merged.

If you wish to read non-traditional file formats (File Per Signal Type or File Per Channel), give the name of the Intan header file (by default, info.rhd or info.rhs) and ensure the corresponding .dat files are present in the same directory.
