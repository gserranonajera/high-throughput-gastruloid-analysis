%% gastruloidAnalysis Example
addpath(genpath('/path/to/utilities')); % Add utility functions to path

microscope_results_path = 'path/to/img_files';
groups = {}; % Leave empty for automatic detection (dish mode)

%% Parameters
param = struct();

% general paramenters
param.plate_mode = false; % true if 96 well plate, false if dishes
param.workers = 20; % number of cores for parallel processing

% imaging channels
param.channel_names = {'BF', 'Bra', 'Sox1'};
param.n_channels = length(param.channel_names);

% parameters for cellpose gastruloid segmentation
param.cellpose_model = 'CP_protrusion_v05'; % pretrained cellpose model to use
param.cellpose_diameter = 140; % diameter for cellpose. Depends on the size of the gastruloids
param.cellpose_downsize = 0.1; % factor for image reduction to segment full gastruloid bodies using cellpose

% parameters for cell probability segmentation
param.probmask.thr = 0.35; % 0.35 minimum probability to find a cell
param.probmask.min_area = 1e4; % minumum area of the mask

%% Tool paths (update for your system)
% image J macro to convert 16bit nd2 images to 8bit tif images
param.macro_8bit_path = 'path/to/macro_convert_single_stack_8bit.ijm';
param.macro_multifile_8bit_path = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\macro_convert_single_stack_8bit_from_multiseries_ND2.ijm';

% python scripts to extract probabilities from cellpose files
param.path2pythonfunction_extract_probs = 'path/to/extract_probabilities_cellpose_seg_object.py';
param.path2pythonfunction_extract_masks = 'path/to/extract_masks_cellpose_seg_object.py';

% python scripts to run morgana
param.morgana_script = 'path/to/run_morgana_pipeline.py';

%% Running the pipeline
gA = gastruloidAnalysis(microscope_results_path, groups, param);

gA.convert8bit();
gA.stepIdx(); gA.save()

gA.focusSurface();
gA.stepIdx(); gA.save()

gA.generateImageChannels();
gA.stepIdx(); gA.save()

gA.findGastruloidBodies();
gA.stepIdx(); gA.save()

gA.findGastruloidCells();
gA.stepIdx(); gA.save()

gA.plotDiagnosis();

gA.extractPropertiesMorgana();
gA.stepIdx(); gA.save()

gA.groupData();
gA.save()

gA.plotSummary()