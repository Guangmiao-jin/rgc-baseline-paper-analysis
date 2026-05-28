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

Raw MEA recordings were collected using the BioCAM DupleX system (3Brain, Lanquart, Switzerland) and saved in `.brw` format.

### Pre-processing output

The Python pre-processing scripts perform spike detection, spike extraction, cluster sorting and stimulation pulse extraction. The final pre-processing output, including spike times and stimulation pulse time trains, is saved in `_cluster.hdf5` format.

### MATLAB analysis

For each recording session, place the corresponding `_cluster.hdf5` file in the single-recording session folder and run `matlabprocessRetinaFlashStimData.m`

This script generates the output structures and plot listed below, as well as PSTH and raster plots for each categorised unit.

- `_responseMetrics.mat` contains extracted response parameters, including bias index, inter-spike interval (ISI) coefficient of variation, ISI violation rate, tau value, peri-stimulus time histogram (PSTH), post-tau firing rate, and the ratio of total spike number during the ON epoch to total spike number during the OFF epoch.

- `_totalneuronsV2.mat`contains the indices of classified cell clusters, including ON transient, ON sustained, OFF transient, OFF sustained, ON–OFF, unconventional and unresponsive units.

- `_psth.mat`contains PSTH values and bin edges for plotting individual cell response traces.
  
- raster plot and PSTH plot for each individual units across 3 light conditions.

Indexing note:
The file names for individual unit plots are based on Python indexing, which starts from 0. However, the indices stored in _totalneuronsV2.mat follow MATLAB indexing, which starts from 1. Therefore, unit numbers in plot file names and MATLAB classification indices differ by 1.

### Plotting functions

After running `processRetinaFlashStimData.m`, the following three `.mat` output files are required for downstream plotting:

- `_responseMetrics.mat`
- `_totalneuronsV2.mat`
- `_psth.mat`

#### Single-recording spatial plots

Run `plot_locations.m` using the corresponding `_cluster.hdf5` and `_totalneuronsV2.mat` files.

This function generates:

- a spatial plot of all classified RGCs projected onto the MEA layout;
- the percentage of responsive units as a function of distance from the centre of the retina, defined here as the optic nerve head;
- a `_distance.mat` file containing the spatial and distance information required for downstream plotting.

#### Optional MEA–retina image alignment

If a live image of the retina on the MEA chip was acquired during recording, the image can be aligned to the MEA grid using `chipGridAlign.m`

This is a semi-manual alignment script. The user selects four reference coordinates on the chip image: top left, top right, bottom left and bottom right. The script then uses the known inter-electrode distance, for example 64 µm in our recording setting, to align the retinal image with the MEA grid.

#### Single-cell PSTH graph under three light conditions

Run `plotAllRasterPSTHs_combined_Response_examples.m` using `plotMetrics.m` which can be derived from `_psth.mat`.

Each `_psth.mat` file contains a `pltcurve` structure with the following fields:

| Field | Type | Description |
|-------|------|-------------|
| `trialPSTHs` | 1×N cell array | Each element corresponds to one neuron, containing a 1×3 cell array representing three light conditions. Each condition stores a matrix of dimensions `(n_repetitions × n_time_bins)`. |
| `PSTH_binEdges` | vector | Time bin edges (in seconds) spanning the pre- and post-stimulation periods. Bin centres are computed as the midpoint between consecutive edges. |

`plotMetrics.m` contains the following fields required:

| Field | Type | Description |
|-------|------|-------------|
| `trialPSTHs` | 1×3 cell array | Each condition stores a matrix of dimensions `(n_repetitions × n_time_bins)`. |
| `PSTH_binEdges` | vector | Time bin edges (in seconds) spanning the pre- and post-stimulation periods. Bin centres are computed as the midpoint between consecutive edges. |

This function generates a PSTH plot with three light conditions for a defined cell

#### Group-level plotting and organising files

For age-group comparisons, place the required output files into folders organised by age group. For example:

```text
rd10_baseline/
├── d23/
│   ├── recording_01_distance.mat
│   ├── recording_01_totalneuronsV2.mat
│   ├── recording_01_psth.mat
│   ├── recording_02_distance.mat
│   ├── recording_02_totalneuronsV2.mat
│   └── recording_02_psth.mat
├── d45/
│   ├── recording_01_distance.mat
│   ├── recording_01_totalneuronsV2.mat
│   └── recording_01_psth.mat
├── d60/
│   ├── recording_01_distance.mat
│   ├── recording_01_totalneuronsV2.mat
│   └── recording_01_psth.mat
...
└── d200/
    ├── recording_01_distance.mat
    ├── recording_01_totalneuronsV2.mat
    └── recording_01_psth.mat
```


After the files are organised by age group, run the following plotting functions:

- `plot_combined_responsiveness.m` plots the overall percentage of responsive units as a function of distance from the retinal centre across all age groups.

- `plot_totalneurons_stacked.m` generates stacked bar plots showing the proportions of classified RGC response types across age groups.

- `psth_onoff_plot.m`plots the median PSTH traces of ON transient, ON sustained, OFF transient and OFF sustained RGCs within a single age group.








