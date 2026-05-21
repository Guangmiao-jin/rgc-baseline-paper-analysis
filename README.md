# Rd10 MEA RGC Response Analysis Pipeline

This repository contains the analysis pipeline used to process rd10 mouse retinal MEA recordings, classify retinal ganglion cell (RGC) light responses, generate summary plots, and compare response subtype percentages across experimental groups.

## Contributors

This analysis pipeline was developed by Guangmiao Jin, with contributions from Dr Michael Savage.

- **Guangmiao Jin** developed the overall pipeline structure, selected and adjusted key analysis and classification parameters, implemented the RGC response classification workflow, and organised the plotting and percentage comparison analyses.
- **Dr Michael Savage** contributed to the Python stimulation pulse extraction script, MATLAB code optimisation, code annotation and commenting, and technical review of the analysis pipeline.

## Repository structure

- **pre-processing**: adapts the HerdingSpikes2 Lightning workflow ([HS2](https://github.com/mhhennig/HS2)), which is based on the SpikeInterface sorting framework. This step performs spike detection, spike extraction, spike sorting and stimulation pulse extraction from raw MEA recordings.
- **analysis**: aligns stimulation time points with spike times, separates noise from electrophysiologically valid sorted units, extracts response metrics, and classifies selected sorted units into ON, OFF, ON–OFF, unconventional and unresponsive response types.
- **plotting**: generates single-recording plots, including sorted unit locations, response subtype percentages, pie charts of classified cell types, and transient/sustained response subtype summaries.

## Input and output data format

### Input data

Raw MEA recordings were collected using the BioCAM DupleX system (3Brain, Lanquart, Switzerland) and saved in `.brw5` format.

### Pre-processing output

The Python pre-processing scripts perform spike detection, spike extraction, cluster sorting and stimulation pulse extraction. The final pre-processing output, including spike times and stimulation pulse time trains, is saved in `.hdf5` format.

### MATLAB analysis

For each recording session, place the corresponding `.hdf5` file in the single-recording session folder and run:

```matlab
processRetinaFlashStimData.m
