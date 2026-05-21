# rd10 MEA RGC Response Analysis Pipeline

This repository contains the analysis pipeline used to process rd10 mouse retinal MEA recordings, classify retinal ganglion cell (RGC) light responses, generate summary plots, and compare response subtype percentages across experimental groups.

## Contributors

This analysis pipeline was developed by Guangmiao Jin, with contributions from Dr Michael Savage.

- **Guangmiao Jin** developed the overall pipeline structure, selected and adjusted key analysis and classification parameters, implemented the RGC response classification workflow, and organised the plotting and percentage comparison analyses.
- **Dr Michael Savage** contributed to the Python stimulation pulse extraction script, MATLAB code optimisation, and technical review of the analysis pipeline.

## Repository structure
- pre-processing : adapt Herding Spikes 2 Lightning (https://github.com/mhhennig/HS2) script which is based on unified sorter SpikeInterface 
- analysis : clean the sorted data by separating the noise from electrophysiologically valid signal, categorize the selected sorted units into ON, OFF, ON-OFF and unconventional types
- plotting: plot the single recording graphs: i.e. the location of the sorted units, the percentage of all cell types in a stacked graph and the proportion of transient/sustained in cells

## Input and output data format

1. Raw data is collected from equipment BioCAM DupleX system (3Brain, Lanquart, Switzerland) and saved in .brw5 format
2. Python script processes spike detection, spike extraction, cluster sorting and stimulation pulse extraction, the final output, including spikes and stimulation pulse time trains is saved in .hdf5 format
3. Put .hdf5 file in a single recording session file, and execute processRetinaFlashStimDataFN1.mat
4. The ouput from the script is: a mat contains all paramters, a mat contains 
