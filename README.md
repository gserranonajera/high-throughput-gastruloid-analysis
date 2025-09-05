# gastruloidAnalysis

**Automated MATLAB pipeline for quantitative analysis of gastruloid and pescoid microscopy data**

A comprehensive workflow for processing time-lapse microscopy images from raw acquisition through statistical analysis and visualization. Designed for developmental biology researchers studying gastruloid development, morphogenesis, and cell migration.

## Features

- **Complete Pipeline**: Raw image conversion → AI segmentation → Feature extraction → Statistical analysis
- **AI-Powered Segmentation**: Cellpose integration for gastruloid body and cell migration detection
- **Advanced Morphological Analysis**: LocoEFA shape analysis, fluorescence profiling, temporal tracking
- **Flexible Formats**: Supports 96-well plates and culture dishes
- **Parallel Processing**: Optimized for large datasets with configurable worker pools
- **Quality Control**: Comprehensive diagnostic visualizations and validation tools

## Workflow Overview

```
Raw .nd2 Files → 8-bit Conversion → Surface Focus Extraction → Channel Generation
       ↓
AI Segmentation (Cellpose) → Feature Extraction (Morgana) → Statistical Analysis
       ↓
PCA Analysis → Time Series → Group Comparisons → Diagnostic Plots
```

## Requirements

### MATLAB Dependencies
- MATLAB R2019b or later
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox  
- Parallel Computing Toolbox (recommended)

### External Software
- **ImageJ/FIJI** - Image preprocessing
- **Python 3.7+** with Cellpose
- **Morgana** - Morphological analysis tools

### Python Environment Setup
```bash
conda create -n cellpose python=3.8
conda activate cellpose
pip install cellpose[gui]

conda create -n morgana python=3.8
conda activate morgana
# Install morgana following their documentation
```

## Installation

1. **Clone Repository**
```bash
git clone https://github.com/username/gastruloidAnalysis.git
cd gastruloidAnalysis
```

2. **Download External Dependencies**
All required utility functions are included in this repository:

**MATLAB Utilities:**
- `openTiffFijiHyperStack.m` - Read multi-dimensional TIFF files
- `openTiffStack.m` - Read TIFF stacks  
- `saveTiffHyperStack.m` - Save multi-dimensional TIFF files
- `saveTiffTimeStack.m` - Save time-series TIFF stacks
- `run_imagej_script.m` - Execute ImageJ macros from MATLAB

**ImageJ Macros:**
- `macro_convert_single_stack_8bit.ijm` - Convert .nd2 to 8-bit TIFF
- `macro_convert_single_stack_8bit_from_multiseries_ND2.ijm` - Extract from multi-series files

**Python Scripts:**
- `extract_probabilities_cellpose_seg_object.py` - Extract probability maps from Cellpose
- `extract_masks_cellpose_seg_object.py` - Extract masks from Cellpose
- `run_morgana_pipeline.py` - Execute Morgana analysis
- `objectsparser.py` - Parse segmentation objects

3. **Configure Paths**
Update the parameter structure in your analysis script or use environment variables:

```matlab
% Option 1: Direct parameter configuration
param = struct();
param.macro_8bit_path = '/path/to/macro_convert_single_stack_8bit.ijm';
param.morgana_script = '/path/to/run_morgana_pipeline.py';
% ... other paths

% Option 2: Environment variables (recommended)
setenv('GASTRULOID_TOOLS_PATH', '/path/to/utility/scripts');
```

## Quick Start

### Basic Analysis Pipeline

```matlab
% Initialize analysis
microscope_path = '/path/to/nd2/files';
groups = {'control', 'treatment1', 'treatment2'};

% Configure parameters (optional)
param = struct();
param.plate_mode = false;  % true for 96-well plates
param.workers = 8;         % parallel workers
param.cellpose_downsize = 0.1;  % image reduction factor

gA = gastruloidAnalysis(microscope_path, groups, param);

% Step 1: Convert to 8-bit
gA.stepIdx();
gA.convert8bit();
gA.save();

% Step 2: Extract surface focus
gA.stepIdx();
gA.focusSurface();
gA.save();

% Step 3: Generate channel projections  
gA.stepIdx();
gA.generateImageChannels();
gA.save();

% Step 4: Segment gastruloids
gA.stepIdx();
gA.findGastruloidBodies();
gA.save();

% Step 5: Detect cell migration (optional)
gA.stepIdx();
gA.findGastruloidCells();
gA.save();

% Step 6: Extract morphological properties
gA.stepIdx();
gA.extractPropertiesMorgana();
gA.save();

% Step 7: Group analysis and visualization
gA.stepIdx();
gA.groupData();
gA.plotSummary();
gA.save();
```

### 96-Well Plate Analysis

```matlab
% For 96-well plate experiments
param.plate_mode = true;
groups = {'01', '02', '03'};  % Column-based grouping

gA = gastruloidAnalysis(microscope_path, groups, param);
% ... run pipeline as above

% Generate plate diagnostic plots
gA.plotDiagnosis96WellPlate(0.1);  % 0.1 = scale factor
```

## Configuration Options

### Key Parameters

```matlab
param = struct();

% Processing settings
param.cellpose_downsize = 0.1;    % Image reduction for segmentation
param.workers = 20;               % Parallel processing cores
param.plate_mode = false;         % 96-well plate format

% Cell probability detection
param.probmask.thr = 0.35;        % Probability threshold
param.probmask.min_area = 1e4;    % Minimum cell area

% File paths (configure for your system)
param.macro_8bit_path = 'path/to/macro_convert_single_stack_8bit.ijm';
param.morgana_script = 'path/to/run_morgana_pipeline.py';
param.path2pythonfunction_extract_probs = 'path/to/extract_probabilities_cellpose_seg_object.py';
param.path2pythonfunction_extract_masks = 'path/to/extract_masks_cellpose_seg_object.py';
```

### File Naming Conventions

**Dish Mode:**
- Format: `001_Treatment1_Sample_001.nd2`
- Groups extracted automatically from "Treatment" patterns

**Plate Mode:**
- Format: `01_A01.nd2` or `A01.nd2`
- Groups based on column numbers or custom patterns

## Advanced Features

### Manual Segmentation Correction
```matlab
% After running findGastruloidBodies(), correct masks in Cellpose GUI
% Then process corrected masks:
gA.findGastruloidManualCorrection();
```

### Custom Morphological Analysis
```matlab
% Analyze specific properties at final timepoint
gA.plotPropertiesBoxPlots(24, {'area', 'spreading_ratio', 'midline_length'});

% PCA analysis of shape evolution
gA.plotPropertiesPCA(-1, true, {'area', 'circularity', 'ch1'});

% LocoEFA shape complexity analysis
gA.plotLocoEfaComplexity(24);
```

### Quality Control
```matlab
% Generate diagnostic montages
gA.plotDiagnosis(0.1);  % For dish experiments
gA.plotDiagnosis96WellPlate(0.1);  % For plate experiments

% Individual sample visualization
gA.plotSingleMask(1, 10);    % Sample 1, timepoint 10
gA.plotSingleImage(1, 10);   % Channel visualization
```

## Output Structure

The pipeline creates organized output directories:

```
analysis_folder/
├── 00_raw/                 # 8-bit converted images
├── 01_surface/            # Surface focus projections  
├── 02_channels/           # Multi-channel projections
├── 03_segmentation/       # Gastruloid body masks
├── 04_cellprob_segmentation/  # Cell migration masks (optional)
├── 05_properties/         # Morgana analysis results
├── 06_group_analysis/     # Statistical analysis and plots
├── mask_diagnosis/        # Quality control images
└── error_log.txt         # Error logging
```

## Data Analysis Features

- **Morphological Properties**: Area, perimeter, eccentricity, axis ratios, form factors
- **Shape Analysis**: LocoEFA coefficients, shape complexity, temporal evolution
- **Fluorescence Quantification**: AP/LR/radial profiles, polarization metrics
- **Cell Migration**: Spreading ratios, migration areas, directionality
- **Statistical Analysis**: PCA, time-series analysis, group comparisons
- **Temporal Tracking**: Consistent labeling across time points

## Troubleshooting

### Common Issues

**Path Errors:**
- Ensure all utility scripts are accessible
- Use absolute paths or properly configure environment variables
- Check that ImageJ and Python environments are properly activated

**Segmentation Quality:**
- Adjust `cellpose_downsize` parameter for your image resolution
- Use diagnostic plots to validate segmentation quality
- Consider manual correction for critical samples

**Memory Issues:**
- Reduce `workers` parameter for large datasets
- Process samples in batches using `indices` parameter
- Monitor MATLAB memory usage during parallel processing

**Platform Compatibility:**
- Windows paths use backslashes, update accordingly
- Conda activation commands may differ between platforms
- Ensure proper file permissions for utility scripts

## Contributing

This research-grade software benefits from community improvements:

1. **Bug Reports**: Use GitHub issues with minimal reproducible examples
2. **Feature Requests**: Describe scientific use cases and expected outputs  
3. **Code Contributions**: Follow MATLAB style guidelines, add unit tests
4. **Documentation**: Improve examples, add troubleshooting guides

## Citation

If you use this software in your research, please cite:

```
Serrano Nájera G, Delahaye A, Steventon B. Gastruloids employ an alternative morphogenetic route to generate a posterior body axis on adherent substrates. bioRxiv. 2025:2025-08.
```

## License

MIT License

MIT License
Copyright (c) 2025 Guillermo Serrano Nájera

## Contact

**Guillermo Serrano Nájera**  
Email: gs714@cam.ac.uk

For technical support:
- Check troubleshooting section
- Review diagnostic output images
- Verify all dependencies are properly installed
- Report issues with complete error logs and system information

---

**Note**: This pipeline was developed for specific microscopy setups and analysis workflows. Adaptation may be required for different experimental conditions, imaging parameters, or analysis requirements.
