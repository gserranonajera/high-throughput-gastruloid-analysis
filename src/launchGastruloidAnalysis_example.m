addpath(genpath('Y:\Users\Guillermo\software\repository'));

microscope_results_path = [pwd filesep 'multipoints'];
groups = {};

cellpose_diameter = 140;
% colormaps
cmap = lines(length(groups));

% dish parameters
param = struct();
param.time_step = 1; % min
param.pixel_sz = 1.3015; % um/pixel
param.cellpose_diameter = cellpose_diameter;
param.cellpose_model = '5X_reduction_plate_dish_v02';
param.cellpose_downsize = 0.2;
param.plate_mode = 0;
param.cmap = cmap;
param.channel_names = {'BF', 'Bra', 'Nuc'};
param.n_channels = length(param.channel_names);
param.link_labels = 0; % add link labels as an option


gA = gastruloidAnalysis(microscope_results_path, groups, param);

gA.convert8bit();
gA.stepIdx(); gA.save()

gA.focusSurface();
gA.stepIdx(); gA.save()

gA.generateImageChannels();
gA.stepIdx; gA.save()

gA.findGastruloidBodies([], [cellpose_diameter], [], [], [], []) % model_name, diameter, indices, time_points, diameter_range, group_name
gA.stepIdx(); gA.save()

gA.findGastruloidCells([], [], [], [], []); %(model_name, diameter, indices, downsize_factor, time_points)
gA.stepIdx(); gA.save()

gA.plotDiagnosis();

gA.extractPropertiesMorgana();
gA.stepIdx(); gA.save()

gA.groupData();
gA.save()

gA.plotSummary()
