# high-throughput-gastruloid-analysis
Automated MATLAB pipeline for gastruloid analysis. AI-powered segmentation, morphological feature extraction, fluorescence quantification, and statistical analysis of time-lapse microscopy data. For developmental biology research.

## Overview

This tool processes raw microscopy images through AI-powered segmentation to extract morphological features, fluorescence patterns, and developmental trajectories. Designed for developmental biology researchers studying gastruloid development and morphogenesis.

## Key Features

- **Automated segmentation** using Cellpose AI models
- **Morphological analysis** including shape descriptors, elongation metrics, and LocoEFA contour analysis from Morgana
- **Fluorescence quantification** with anterior-posterior and radial profiling  
- **Statistical analysis** with PCA, time-series plotting, and group comparisons
- **Multi-format support** for both 96-well plates and culture dishes
- **Quality control** with comprehensive diagnostic visualizations

## Workflow

```
Raw Images → Preprocessing → AI Segmentation → Feature Extraction → Statistical Analysis → Visualization
```

## Requirements

### MATLAB Dependencies
- MATLAB R2019b or later
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (recommended for performance)

### External Software Dependencies

#### ImageJ/FIJI
- **ImageJ/FIJI** - Required for image preprocessing
- **Custom ImageJ Macros** (not included):
  - `macro_convert_single_stack_8bit.ijm`
  - `macro_convert_single_stack_8bit_from_multiseries_ND2.ijm` 
  - `macro_stitch_CNSWE.ijm`

#### Python Environment
- **Python 3.7+**
- **Cellpose** - AI segmentation models
```bash
conda create -n cellpose python=3.8
conda activate cellpose
pip install cellpose[gui]
```

#### External Analysis Tools
- **Morgana analysis functions** (external dependency)
- After installing Morgana the ImageTools/objectsparser.py needs to be substituted by the objectsparser.py provided in this repository to deal with more than 1 gastruloid per image.
  - `run_morgana_analysis_flourescence()`
  - `focus_Nikon_experiments_multistack()`
  - `run_gastruloid_body_cell_detection()`


## Installation

1. **Clone this repository**
```bash
git clone https://github.com/username/gastruloidAnalysis.git
cd gastruloidAnalysis
```

2. **Set up Python environment**
```bash
conda create -n cellpose python=3.8
conda activate cellpose
pip install cellpose[gui]
```

3. **Install ImageJ/FIJI** and required macros

4. **Update file paths** in `gastruloidAnalysis.m`:
   - Modify hard-coded paths to match your system
   - Update macro paths in constructor

## Quick Start

```matlab
% Initialize analysis
microscope_path = 'path/to/your/data';
groups = {'control', 'treatment1', 'treatment2'};
gA = gastruloidAnalysis(microscope_path, groups);

% Run preprocessing
gA.stepIdx(); % Step 1: Raw processing
gA.convert8bit();
gA.stepIdx(); gA.save();

% Surface focus extraction  
gA.focusSurface();
gA.stepIdx(); gA.save();

% Generate channels
gA.generateImageChannels();
gA.stepIdx(); gA.save();

% Segmentation
gA.findGastruloids();
gA.stepIdx(); gA.save();

% Analysis
gA.extractPropertiesMorgana();
gA.stepIdx(); gA.save();

% Group analysis and visualization
gA.groupData();
gA.plotSummary();
```

## Configuration

Key parameters can be modified in the constructor:
```matlab
param.time_step = 60; % minutes between frames
param.pixel_sz = 1.3015; % μm/pixel
param.cellpose_diameter = 30; % pixels (30 for D4, 15 for D2)
param.cellpose_model = '10xreductionmodel';
param.cellpose_downsize = 0.1; % reduction factor
param.plate_mode = false; % true for 96-well plates
```

## File Structure

Input files should follow naming conventions:
- **Dish mode**: `001_Treatment1_Sample_001.tif`
- **Plate mode**: `01_A01.tif` (well format)

## Output

The pipeline generates:
- Segmentation masks
- Morphological measurements
- Fluorescence profiles
- Statistical plots (PCA, box plots, time series)
- Quality control visualizations

## Limitations

- **Platform dependency**: Hard-coded Windows paths need modification
- **External dependencies**: Requires multiple software installations
- **Research-grade software**: May need adaptation for different experimental setups
