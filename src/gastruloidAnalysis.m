classdef gastruloidAnalysis < handle
    %}
    %}

    properties
        main_dir
        microscope_results_path
        raw_dir
        surf_dir
        channel_dir
        segm_dir
        cellprob_segm_dir
        semant_dir % deprecated
        properties_dir % properties from MorGana
        group_analysis_dir %result folder
        examples_dir %not used at the moment

        sample_names
        groups % numerical indices. Experimental groups for plotting and statistics
        group_final_names % Actual string names
        

        error_path
        workers % number of workers for parallel processing
        current_idx % index of the current analysis folder

        macro_8bit_path % image j macro
        macro_stitch_path % image j macro

        group_analysis_param % parameters for group analysis
        param % parameters for analysis
    end

    methods

        % Constructor that initializes the analysis pipeline with data paths, experimental groups, and analysis parameters
        function gA = gastruloidAnalysis(microscope_results_path, groups, param, main_dir)
            
            gA.microscope_results_path = microscope_results_path;

            if ~exist('main_dir','var') || isempty(main_dir)
                gA.main_dir=pwd;
            else
                gA.main_dir=[main_dir filesep];
            end

            % default parameters
            if ~exist('param','var') || isempty(param)
                gA.param = struct();
                gA.param.time_step = 60; % min
                gA.param.pixel_sz = 1.3015; % um/pixel
                gA.param.cellpose_diameter = 30; % (30 D4) (15 D2)
                gA.param.cellpose_model = '10xreductionmodel';
                gA.param.cellpose_downsize = 0.1; % (0.1 D4) (0.2 D2)
                gA.param.plate_mode = false; % true if 96 well plate, false if dishes
            else
                gA.param=param;
            end

            gA.current_idx = 0; % normally idx = 0 is raw
            gA.workers = 20;
            gA.error_path = [gA.main_dir filesep 'error_log.txt'];
    
            gA.getSampleNames(); % populates the property sample_names

            if ~exist('groups','var') || isempty(groups)
                gA.extractGroupNames();
                if isempty(gA.param.cmap)
                    gA.param.cmap = lines(length(gA.groups));
                end
            else
                gA.groups = groups;
            end
            
            % in plate mode samples have names as in a 96well plate A01, A02...
            if gA.param.plate_mode == false; gA.extractGroupFinalNames; end
            
            % image J macro to convert 16bit nd2 images to 8bit tif images
            gA.macro_8bit_path = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\macro_convert_single_stack_8bit.ijm';
            
            % deprecated
            gA.group_analysis_param = struct();
            gA.group_analysis_param.time_step = gA.param.time_step;
            gA.group_analysis_param.max_dist_centre = 1000;
            gA.group_analysis_param.min_area = 300;
        end

        function stepIdx(gA)
            gA.current_idx = gA.current_idx + 1;
        end

        function save(gA)
            save('gA')
        end

        
        %% Initialisation functions for treatment names and indices
        % Extracts sample names from .tif or .nd2 files in the raw directory
        function getSampleNames(gA)

            raw_stacks = dir([gA.raw_dir filesep '*.tif']);
            if length(raw_stacks) == 0
                raw_stacks = dir([gA.microscope_results_path filesep '*.nd2']);
            end
            gA.sample_names = cell(1, length(raw_stacks));

            for f_id = 1:length(raw_stacks)
                name = strsplit(raw_stacks(f_id).name, '.');
                name = name{1};
                gA.sample_names{f_id} = name;
            end
        end

        % Parses treatment condition names from filename using regex patterns
        function extracted_string = extractTreatmentName(gA, input_string)
            % Define the regular expression patterns
            
            treatment_pattern = '\d+_Treatment\d+_(\w+)_\d+';
            % treatment_pattern = '\d+_Treatment\d+_(\w+)';
            % treatment_pattern = '\d+_(*)_\d+';
            control_pattern = '\d+_(\w+)_\d+';
            
            % Try to match the treatment pattern first
            treatment_match = regexp(input_string, treatment_pattern, 'tokens', 'once');
            
            if ~isempty(treatment_match)
                % If treatment pattern matches, return the extracted string
                extracted_string = treatment_match{1};
            else
                % If treatment pattern doesn't match, try the control pattern
                control_match = regexp(input_string, control_pattern, 'tokens', 'once');
                
                if ~isempty(control_match)
                    
                    % If control pattern matches, return the extracted string
                    extracted_string = control_match{1};
                else
                    % If neither pattern matches, return an empty string
                    extracted_string = '';
                    
                end
            end
        end
        
        % Extracts treatment numbers or "control" labels from filenames
        function extracted_string = extractTreatmentNumber(gA, input_string)
            % Define the regular expression pattern
            pattern = '(Treatment\d+)';

            % Use regexp to find matches
            match = regexpi(input_string, pattern, 'tokens', 'once');
            
            % Extract the string if a match is found
            if ~isempty(match)
                extracted_string = match{1};
            else
                pattern = '(control)';
                match = regexpi(input_string, pattern, 'tokens', 'once');
                extracted_string = match{1};
            end
        end

        % Automatically determines experimental groups from sample filenames
        function extractGroupNames(gA)
            gA.groups = unique(cellfun(@(x) gA.extractTreatmentNumber(x), gA.sample_names, 'UniformOutput', false), 'stable');
            % second_elem = @(x) x{2};
            % gA.group_final_names = unique(cellfun(@(x) second_elem(strsplit(x, '_')), gA.sample_names, 'UniformOutput', false), 'stable');
        end
        
        % Creates final group names for display purposes (dish mode only)
        function extractGroupFinalNames(gA)
            gA.group_final_names = unique(cellfun(@(x) gA.extractTreatmentName(x), gA.sample_names, 'UniformOutput', false), 'stable');
            % second_elem = @(x) x{2};
            % gA.group_final_names = unique(cellfun(@(x) second_elem(strsplit(x, '_')), gA.sample_names, 'UniformOutput', false), 'stable');

        end
        
        %% Image Conversion function
        % Combines .nd2 files from two folders with reindexed names. For
        % example when you image in 2 batches and one goes from 1:30 and
        % the other from 1:40, but you need 1:70.
        function mergeMultifolderExperiment(gA, folder_A, folder_B)
            get_first = @(x) x{1};

            A_data = dir([folder_A filesep '*.nd2']);
            A_indices = cell2mat(cellfun(@(x) str2double(get_first(strsplit(x,'_'))), {A_data.name}, 'UniformOutput', false));
            last_A_idx = max(A_indices);


            B_data = dir([folder_B filesep '*.nd2']);
            for idx = 1:length(B_data)
                source = [B_data(idx).folder filesep B_data(idx).name];

                split_name = strsplit(B_data(idx).name, '_');
                newidx = num2str(str2double(get_first(split_name))+last_A_idx);
                modified_idx_name = [newidx '_' split_name{2} '_' split_name{3}];

                detination = [folder_A filesep modified_idx_name];

                copyfile(source, detination)
            end
            
%             B_indices = cellfun(@(x) str2double(get_first(strsplit(x,'_'))), {B_data.name}, 'UniformOutput', false);
        end

        % Converts raw microscopy files in the "multipoints" folder to 8-bit format using ImageJ macros
        function convert8bit(gA, indices)
            
            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            suffix = '.nd2';
            output_folder = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_raw'];
            mkdir(output_folder);

            % initialise the parallel pool
            pool = gcp('nocreate');
            if isempty(pool)
               parpool(gA.workers);
            elseif pool.NumWorkers~=gA.workers
               delete(pool);
               parpool(gA.workers);
            end

            parfor s_id = indices
                input_file_path = [gA.microscope_results_path filesep  gA.sample_names{s_id} suffix];
                run_imagej_script(gA.macro_8bit_path, [], [], input_file_path, output_folder, suffix);
            end

            gA.raw_dir = output_folder;
        end

        % Extracts and converts individual series from multi-series
        % "multipoints.nd2"
        % files. Use when individual images in "multipoints" folder are corrupted.
        function convert8bitFromMultiseriesFile(gA, indices, big_file_path)
            
            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            output_folder = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_raw'];
            mkdir(output_folder);

            gA.macro_8bit_path = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\macro_convert_single_stack_8bit_from_multiseries_ND2.ijm';

            pool = gcp('nocreate');
            if isempty(pool)
               parpool(gA.workers);
            elseif pool.NumWorkers~=gA.workers
               delete(pool);
               parpool(gA.workers);
            end

            parfor s_id = indices
                img_name = gA.sample_names{s_id};
                run_imagej_script(gA.macro_8bit_path, [], [], big_file_path, output_folder, num2str(s_id), img_name);
            end

            gA.raw_dir = output_folder;
        end

        %% Stitch functions for broken stage
        % used when stage was broken
        function stitchImagesMacro(gA, indices)           

            raw_stacks = dir([gA.raw_dir filesep '*_C.tif']);
%             raw_stacks = dir([gA.raw_dir filesep '*.tif']);
            gA.sample_names = cell(1, length(raw_stacks));

            for f_id = 1:length(raw_stacks)
                name = extractBefore(raw_stacks(f_id).name, '_C.tif');
                %name = extractBetween(raw_stacks(f_id).name, digitsPattern(3), '_C.tif');
                gA.sample_names{f_id} = name;
            end

            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end
            
            output_folder = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_stitch'];
            mkdir(output_folder);

            gA.macro_stitch_path = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\macro_stitch_CNSWE.ijm';

            pool = gcp('nocreate');
            if isempty(pool)
               parpool(gA.workers);
            elseif pool.NumWorkers~=gA.workers
               delete(pool);
               parpool(gA.workers);
            end

            parfor s_id = indices
                img_name = extractAfter(gA.sample_names{s_id}, digitsPattern(3)+'_');
                index = str2double(extractBefore(gA.sample_names{s_id}, '_'));
                name_C  = [num2str(index+0, '%03d') '_' img_name '_C'];
                name_NE = [num2str(index+1, '%03d') '_' img_name '_NE'];
                name_NW = [num2str(index+2, '%03d') '_' img_name '_NW'];
                name_SW = [num2str(index+3, '%03d') '_' img_name '_SW'];
                name_SE = [num2str(index+4, '%03d') '_' img_name '_SE'];
                run_imagej_script(gA.macro_stitch_path, [], [],...
                                  gA.raw_dir,... 
                                  name_C,...
                                  name_NE,...
                                  name_NW,...
                                  name_SW,...
                                  name_SE,...
                                  output_folder,...
                                  img_name);
            end

            gA.raw_dir = output_folder;
        end

        % used when stage was broken
        function simpleStitch(gA)

            raw_stacks = dir([gA.raw_dir filesep '*_C.tif']);
            gA.sample_names = cell(1, length(raw_stacks));
            stitch_folder = [gA.main_dir filesep num2str(gA.current_idx, '%02d')  '_raw_stitch'];
            mkdir(stitch_folder)

            for f_id = 1:length(raw_stacks)
                name = extractBetween(raw_stacks(f_id).name, digitsPattern(3)+'_', '_C.tif');
                gA.sample_names{f_id} = name{1};
            end

            parfor s_idx = 1:length(gA.sample_names)
                s_idx
                tiles = {dir([gA.raw_dir filesep '*' gA.sample_names{s_idx} '*']).name};
                C  = openTiffFijiHyperStack([gA.raw_dir filesep tiles{1}]);

                NE = zeros(size(C), 'uint8');
                NW = zeros(size(C), 'uint8');
                SW = zeros(size(C), 'uint8');
                SE = zeros(size(C), 'uint8');

                try
                    NE = openTiffFijiHyperStack([gA.raw_dir filesep tiles{2}]);
                    NW = openTiffFijiHyperStack([gA.raw_dir filesep tiles{3}]);
                    SW = openTiffFijiHyperStack([gA.raw_dir filesep tiles{4}]);
                    SE = openTiffFijiHyperStack([gA.raw_dir filesep tiles{5}]);
                end

                vol = zeros(size(C, 1)*2, size(C, 2)*2, size(C, 3), size(C, 4), size(C, 5), 'uint8');
                vol(1:size(C), 1:size(C), :, :, :) = NE;
                vol(1:size(C), size(C)+1:end, :, :, :) = NW;
                vol(size(C)+1:end, size(C)+1:end, :, :, :) = SW;
                vol(size(C)+1:end,1:size(C), :, :, :) = SE;
                vol(size(C)/2:size(C)*3/2-1, size(C)/2:size(C)*3/2-1, :, :, :) = C;

                %dimensions: [x, y, z, t, c]
                fiji_descr = ['ImageJ=1.52p' newline ...
                            'images=' num2str(size(vol,3)*size(vol,4)*size(vol,5)) newline... 
                            'slices=' num2str(size(vol,3)) newline...
                            'frames=' num2str(size(vol,4)) newline...
                            'channels=' num2str(size(vol,5)) newline...
                            'hyperstack=true' newline...
                            'mode=8bit' newline...  
                            'loop=false' newline...  
                            'min=0.0' newline...      
                            'max=255'];  % 8bit image
                
                t = Tiff([stitch_folder filesep num2str(s_idx, '%03d') '_' gA.sample_names{s_idx} '.tif'], 'w');
                tagstruct = {};
                tagstruct.ImageLength = size(vol,1);
                tagstruct.ImageWidth = size(vol,2);
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;%Tiff.Photometric.RGB; % rgb
                tagstruct.BitsPerSample = 8; % 8bit
                tagstruct.SamplesPerPixel = 1;% 3rgb
                tagstruct.Compression = Tiff.Compression.LZW;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                tagstruct.ImageDescription = fiji_descr;

                % frame, slice, channel order is important
                for frame = 1:size(vol,4)
                    for slice = 1:size(vol,3)
                        for channel = 1:size(vol,5)
                                
                            t.setTag(tagstruct)
                            t.write((vol(:,:,slice,frame,channel)));
                            t.writeDirectory(); % saves a new page in the tiff file
                        end
                    end
                end
                t.close() 
            end

            gA.raw_dir = stitch_folder;

        end
        
        % used when stage was broken
        function cropImages(gA, cube_limits, indices)
            %cube_limits = [xmin ymin zmin width height depth]
            xmin = cube_limits(1);
            ymin = cube_limits(2);
            zmin = cube_limits(3);
            width = cube_limits(4);
            height = cube_limits(5);
            depth = cube_limits(6);
            
            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end
            
            parfor s_id = indices
                try
                    filename = [gA.raw_dir filesep gA.sample_names{s_id} '.tif'];
                    img_vol_full = openTiffFijiHyperStack(filename); % x, y, x, t, c
                    vol = img_vol_full(xmin:xmin+width-1, ymin:ymin+height-1, zmin:zmin+depth-1, :, :);

                    %dimensions: [x, y, z, t, c]
                    fiji_descr = ['ImageJ=1.52p' newline ...
                                'images=' num2str(size(vol,3)*size(vol,4)*size(vol,5)) newline... 
                                'slices=' num2str(size(vol,3)) newline...
                                'frames=' num2str(size(vol,4)) newline...
                                'channels=' num2str(size(vol,5)) newline...
                                'hyperstack=true' newline...
                                'mode=8bit' newline...  
                                'loop=false' newline...  
                                'min=0.0' newline...      
                                'max=255'];  % 8bit image
                    
                    t = Tiff(filename, 'w');
                    tagstruct = {};
                    tagstruct.ImageLength = size(vol,1);
                    tagstruct.ImageWidth = size(vol,2);
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;%Tiff.Photometric.RGB; % rgb
                    tagstruct.BitsPerSample = 8; % 8bit
                    tagstruct.SamplesPerPixel = 1;% 3rgb
                    tagstruct.Compression = Tiff.Compression.LZW;
                    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                    tagstruct.ImageDescription = fiji_descr;

                    % frame, slice, channel order is important
                    for frame = 1:size(vol,4)
                        for slice = 1:size(vol,3)
                            for channel = 1:size(vol,5)
                                    
                                t.setTag(tagstruct)
                                t.write((vol(:,:,slice,frame,channel)));
                                t.writeDirectory(); % saves a new page in the tiff file
                            end
                        end
                    end
                    t.close() 

                catch ME
                    fprintf(1,'The identifier was:\n%s', ME.identifier);
                    fprintf(1,'There was an error! The message was:\n%s', ME.message);
      
                    fid = fopen(gA.error_path,'a+');
                    fprintf(fid,'%s\n',ME.message);
                    for e=1:length(ME.stack)
                      fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                    end
                    
                    % close file
                    fclose(fid);
                end
            end

        end


        %% pre-processing functions
        % Extracts focused surface from 3D z-stacks
        function focusSurface(gA, indices, shift_idx)

            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            if ~exist('shift_idx', 'var') || isempty(shift_idx)
                shift_idx = 0;
            end

            gA.surf_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_surface'];
            mkdir(gA.surf_dir);

            for s_id = indices
                try
                    filename = [gA.raw_dir filesep gA.sample_names{s_id} '.tif'];
                    outdir = [gA.surf_dir filesep ];
                    focus_Nikon_experiments_multistack(filename, outdir, gA.sample_names{s_id}, gA.workers, shift_idx)
                catch ME
                    fprintf(1,'The identifier was:\n%s', ME.identifier);
                    fprintf(1,'There was an error! The message was:\n%s', ME.message);
      
                    fid = fopen(gA.error_path,'a+');
                    fprintf(fid,'%s\n',ME.message);
                    for e=1:length(ME.stack)
                      fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                    end
                    
                    % close file
                    fclose(fid);
                end
            end
        end
        
        % Creates 2D projections combining surface data with channel max projections
        function generateImageChannels(gA, indices)

            if ~exist('indices','var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            gA.channel_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_channels'];
            mkdir(gA.channel_dir);

            parfor s_id = indices
                try

                    S = openTiffFijiHyperStack([gA.raw_dir filesep gA.sample_names{s_id} '.tif']);
                    I = openTiffFijiHyperStack([gA.surf_dir filesep gA.sample_names{s_id} '.tif']);

                    new_stack = zeros(size(S,1), size(S,2), 1, size(S,4), size(S,5), 'uint8');
                    for t_idx = 1:size(S, 4)

                        new_stack(:,:,1, t_idx,1) = I(:,:,1,t_idx); % surface first
                        for ch_idx = 2:size(S, 5) % add other channels
                            new_stack(:, :, 1, t_idx, ch_idx) =uint8(max(squeeze(S(:, :, :, t_idx, ch_idx)),[],3)); % maximum projection;
                        end
                    end

                    outpath = [gA.channel_dir filesep gA.sample_names{s_id} '.tif'];
                    saveTiffHyperStack(new_stack, outpath)

                catch ME
                    fprintf(1,'The identifier was:\n%s', ME.identifier);
                    fprintf(1,'There was an error! The message was:\n%s', ME.message);
      
                    fid = fopen(gA.error_path,'a+');
                    fprintf(fid,'%s\n',ME.message);
                    for e=1:length(ME.stack)
                      fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                    end
                    
                    % close file
                    fclose(fid);
                end
            end

        end

        
        %% Segmentation functions
        % Segments gastruloids using Cellpose AI with specified model and
        % diameter (deprecated)
        function findGastruloids(gA, model_name, diameter, indices) 

            if ~exist('model_name', 'var') || isempty(model_name)
                model_name = gA.param.cellpose_model;
            end
            
            if ~exist('diameter', 'var') || isempty(diameter)
                diameter = gA.param.cellpose_diameter;
            end

            if ~exist('indices', 'var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            downsize_factor = gA.param.cellpose_downsize;
            gA.segm_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_segmentation'];
            reduced_dir = [gA.segm_dir filesep 'reduced_images'];
            raw_mask_dir = [gA.segm_dir filesep 'raw_mask'];
            mkdir(gA.segm_dir);
            mkdir(reduced_dir)

            parfor s_idx_idx = 1:length(indices)
                s_idx = indices(s_idx_idx);
                disp(s_idx)
                S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);

                for t_idx = 1:size(S1, 3)
                    imwrite(uint8(imresize(S1(:,:,t_idx), downsize_factor)), [reduced_dir filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '.tif'])
                end
                
%                 S2 = zeros([ceil(size(S1,1)*downsize_factor), ceil(size(S1,2)*downsize_factor), size(S1,3)]);
%                 
%                 for z = 1:size(S1, 3)
%                     S2(:,:,z) =  imresize(S1(:,:,z), downsize_factor);
%                 end
                

%                 saveTiffStack(uint8(S2), [reduced_dir filesep gA.sample_names{s_idx} '.tif'])
            end


            try

                 command = ['conda activate cellpose &',...
                    'python -m cellpose --pretrained_model ' model_name, ...
                    ' --diameter ' num2str(diameter), ...
                    ' --dir ' reduced_dir, ...
                    ' --savedir ' raw_mask_dir, ...
                    ' --use_gpu', ...
                    ' --save_tif',...
                    ' --stitch_threshold 0.0',...
                    ' --fast_mode',...
                    ];
                disp(command)
                system(command)


                % this shoud be moved out of this loop
                parfor s_idx_idx = 1:length(indices) % could be in paralle
                    s_idx = indices(s_idx_idx);

                    S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                    final_dim = size(S1);
                    
                    path2label = [raw_mask_dir filesep gA.sample_names{s_idx} '*_cp_masks.tif'];
                    linkLabels(path2label, gA.segm_dir, [gA.sample_names{s_idx} '_mask.tif'], [], final_dim, gA.param.link_labels);
%                     delete(path2label)
                end

            catch ME
                fprintf(1,'The identifier was:\n%s', ME.identifier);
                fprintf(1,'There was an error! The message was:\n%s', ME.message);
  
                fid = fopen(gA.error_path,'a+');
                fprintf(fid,'%s\n',ME.message);
                for e=1:length(ME.stack)
                  fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                end
                
                % close file
                fclose(fid);
            end

        end

        % Time-specific gastruloid body segmentation with variable diameter support
        function findGastruloidBodies(gA, model_name, diameter, indices, time_points, diameter_range, group_name)

            if ~exist('model_name', 'var') || isempty(model_name)
                %model_name = 'only_body_nospreading_v04';
                %model_name = 'CP_test_protrusion_v04';
                model_name = 'CP_protrusion_v05';
            end
            
            if ~exist('diameter', 'var') || isempty(diameter)
                diameter = 200;
            end

            if ~exist('indices', 'var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            if ~exist('diameter_range', 'var') || isempty(diameter_range)
               diameter_range = [diameter, diameter];
            end

            if ~exist('time_points', 'var') || isempty(time_points)
                S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
               time_points = 1:size(S, 3);
            end

            if ~exist('group_name', 'var') || isempty(group_name)
                group_name = [];
            end

            diameter_vector = linspace(diameter_range(1), diameter_range(2), length(time_points));

            downsize_factor = gA.param.cellpose_downsize;
            gA.segm_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_segmentation'];
            raw_mask_dir = [gA.segm_dir filesep 'raw_mask'];
            mkdir(gA.segm_dir);      
            
            
            for t_idx_idx = 1:length(time_points)
                t_idx = time_points(t_idx_idx);
                reduced_dir = [gA.segm_dir filesep 'reduced_images' filesep group_name filesep num2str(t_idx)];
                mkdir(reduced_dir)

                parfor s_idx_idx = 1:length(indices)
                    s_idx = indices(s_idx_idx);
                    disp(s_idx)
                    S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                    imwrite(uint8(imresize(S1(:,:,t_idx), downsize_factor)), [reduced_dir filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '.tif'])
                    
                end


                if ~isempty(diameter_range)

                    diameter = diameter_vector(t_idx_idx);
                end

                try
    
                     command = ['conda activate cellpose &',...
                        ' python -m cellpose --pretrained_model ' model_name, ...
                        ' --diameter ' num2str(diameter), ...
                        ' --dir ' reduced_dir, ...
                        ' --savedir ' raw_mask_dir, ...
                        ' --use_gpu', ...
                        ' --save_tif',...
                        ' --stitch_threshold 0.0',...
                        ' --fast_mode',...
                        ];
    
                    disp(command)
                    system(command)

%                     rmdir(reduced_dir, 's')
    
                catch ME
                    fprintf(1,'The identifier was:\n%s', ME.identifier);
                    fprintf(1,'There was an error! The message was:\n%s', ME.message);
      
                    fid = fopen(gA.error_path,'a+');
                    fprintf(fid,'%s\n',ME.message);
                    for e=1:length(ME.stack)
                      fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                    end
                    
                    % close file
                    fclose(fid);
                end

            end % timepoints

            parfor s_idx_idx = 1:length(indices) % this feels wrong. this should be done only once after the time loop CHECK FUTURE GUILLERMO
                s_idx = indices(s_idx_idx)

                S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                final_dim = size(S1);
                
                path2label = [raw_mask_dir filesep gA.sample_names{s_idx} '*_cp_masks.tif'];
                linkLabels(path2label, gA.segm_dir, [gA.sample_names{s_idx} '_mask.tif'], [], final_dim, gA.param.link_labels);
%                         delete(path2label)
            end

        end

        % Generates cell probability maps using Cellpose for migration analysis
        function findGastruloidCells(gA, model_name, diameter, indices, downsize_factor, time_points)

            %{
            Using CP model with full resolution 20X image to find the
            probability maps for individual cells. CP model has been
            pre-trained to segment cells and the probability maps it
            produces gives a good approximation of the the areas where the
            cells are migrating out of the gastruloid.

            For this reason the model_name and diameter should be mainted
            as default. It does not make sense to change the diameter,
            because the radii of cells does not change, unlike the diameter
            of the gastruloids.
            %}

            if ~exist('model_name', 'var') || isempty(model_name)
                model_name = 'CP';
            end
            
            if ~exist('diameter', 'var') || isempty(diameter)
                diameter = 15;
            end

            if ~exist('indices', 'var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            if ~exist('downsize_factor', 'var') || isempty(downsize_factor)
                downsize_factor = 1; % no resize
            end
            
            if ~exist('time_points', 'var') || isempty(time_points) 
                S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
                time_points = [1:size(S, 3)];
            end

            
            gA.cellprob_segm_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_cellprob_segmentation'];
            reduced_dir = [gA.cellprob_segm_dir filesep 'reduced_images'];
            raw_mask_dir = [gA.cellprob_segm_dir filesep 'raw_mask'];
            mkdir(gA.cellprob_segm_dir);
            mkdir(reduced_dir)
            mkdir(raw_mask_dir)

            parfor s_idx_idx = 1:length(indices)
                s_idx = indices(s_idx_idx);
                disp(s_idx)
                try
                    S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
    
                    for t_idx = time_points
                        imwrite(uint8(imresize(S1(:,:,t_idx), downsize_factor)), [reduced_dir filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '.tif'])
                    end
                
                catch
                    continue
                end
                
            end

            try

                 command = ['conda activate cellpose &',...
                    'python -m cellpose --pretrained_model ' model_name, ...
                    ' --diameter ' num2str(diameter), ...
                    ' --dir ' reduced_dir, ...
                    ' --savedir ' raw_mask_dir, ...
                    ' --use_gpu', ...
                    ' --stitch_threshold 0.0',...
                    ' --fast_mode',...
                    ];
                disp(command)
                system(command)
                

                gA.generateCellProbMask()

                parfor s_idx_idx = 1:length(indices) % could be in paralle
                    s_idx = indices(s_idx_idx);

                    S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                    final_dim = size(S1);
                    
                    path2label = [raw_mask_dir filesep gA.sample_names{s_idx} '*mask.tif'];
                    linkLabels(path2label, gA.cellprob_segm_dir, [gA.sample_names{s_idx} '_mask.tif'], [], final_dim, gA.param.link_labels);
%                     delete(path2label)
                end
            
            catch ME
                fprintf(1,'The identifier was:\n%s', ME.identifier);
                fprintf(1,'There was an error! The message was:\n%s', ME.message);
  
                fid = fopen(gA.error_path,'a+');
                fprintf(fid,'%s\n',ME.message);
                for e=1:length(ME.stack)
                  fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                end
                
                % close file
                fclose(fid);
            end

        end

        % Extracts probability maps from Cellpose segmentation objects using a
        % python function since .npy cellpose files can only be opened that
        % way
        function extractProbabilitiesCellpose(gA, infolder, outfolder)

            path2pythonfunction = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\extract_probabilities_cellpose_seg_object.py';
            filenames = dir([infolder filesep '*_seg.npy']);

            parfor s_idx = 1:length(filenames)
                filename = strsplit(filenames(s_idx).name, "_seg.npy");
                filename = filename{1};
                command = [ 'conda activate cellpose & ',...
                            'python ' path2pythonfunction, ...
                            ' ' filename,...
                            ' ' infolder,...
                            ' ' outfolder,...
                            ];
                disp(command)
                system(command)
            end

        end

        % Converts cell probability maps to binary masks using thresholding
        function generateCellProbMask(gA)
            infolder = [gA.cellprob_segm_dir filesep 'reduced_images'];
            outfolder = [gA.cellprob_segm_dir filesep 'probabilities'];
            mkdir(outfolder)
            gA.extractProbabilitiesCellpose(infolder, outfolder)
            
            %%
            
            prob_imgs = dir([outfolder filesep '*.tif']);
            
            thr = 0.35; %0.35 minimum prob. to find a cell
            min_area = 1e4;
            
            parfor s_idx = 1:length(prob_imgs)
                filename = [prob_imgs(s_idx).folder filesep prob_imgs(s_idx).name];
                base_name = strsplit(prob_imgs(s_idx).name, 'prob_');
                base_name = base_name{2};
                
                I = imread(filename);
                P = double(I)./255;
                
                % imagesc(P);
                % axis image
                % colorbar
                
                BW = P>thr; % binary image
                BW = bwareaopen(BW, min_area);
            %     BW = imclose(BW, strel('disk', 15));
            %     BW = imgaussfilt(double(BW), 5) > 0.5;
                BW = imfill(BW, 'holes');
                % imshow(BW)
                
                bound = imdilate(bwperim(BW), strel('disk', 5));
                IB = imread([infolder filesep base_name]);
                % imshow(imfuse(imadjust(IB), bound))
            %     title(base_name)
                imwrite(BW, [gA.cellprob_segm_dir filesep 'raw_mask' filesep base_name '_mask.tif'])
            end
        end

        function generateBodyCellMask(gA) % deprecated

            gA.semant_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_semantic_segmentation'];
            mkdir(gA.semant_dir);
            
            try
                run_gastruloid_body_cell_detection(gA.surf_dir, gA.semant_dir)
            
            catch ME
                fprintf(1,'The identifier was:\n%s', ME.identifier);
                fprintf(1,'There was an error! The message was:\n%s', ME.message);
  
                fid = fopen(gA.error_path,'a+');
                fprintf(fid,'%s\n',ME.message);
                for e=1:length(ME.stack)
                  fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                end
                
                % close file
                fclose(fid);
            end

        end

        %% Manual correction functions
        function extractMaskCellpose(gA, infolder, outfolder)
            % intended for manual corrections

            path2pythonfunction = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\extract_masks_cellpose_seg_object.py';
            filenames = dir([infolder filesep '*_seg.npy']);

            parfor s_idx = 1:length(filenames)
                filename = strsplit(filenames(s_idx).name, "_seg.npy");
                filename = filename{1};
                command = [ 'conda activate cellpose & ',...
                            'python ' path2pythonfunction, ...
                            ' ' filename,...
                            ' ' infolder,...
                            ' ' outfolder,...
                            ];
                disp(command)
                system(command)
            end

        end

        % Processes manually corrected Cellpose segmentations
        function findGastruloidManualCorrection(gA, indices, time_points, group_name)
            % correct with cellpose GUI and then run this.
            
            if ~exist('indices', 'var') || isempty(indices)
                indices = 1:length(gA.sample_names);
            end

            if ~exist('time_points', 'var') || isempty(time_points)
                time_points = [];
            end

            if ~exist('group_name', 'var') || isempty(group_name)
                group_name = [];
            end


            raw_mask_dir = [gA.segm_dir filesep 'raw_mask'];

            if isempty(time_points)
                reduced_dir = [gA.segm_dir filesep 'reduced_images'];
                gA.extractMaskCellpose(reduced_dir, raw_mask_dir)
    
                try
    
                    parfor s_idx_idx = 1:length(indices)
                        s_idx = indices(s_idx_idx)
    
                        S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                        final_dim = size(S1);
                        
                        path2label = [raw_mask_dir filesep gA.sample_names{s_idx} '*_cp_masks.tif'];
                        linkLabels(path2label, gA.segm_dir, [gA.sample_names{s_idx} '_mask.tif'], [], final_dim, gA.param.link_labels);
    %                         delete(path2label)
                    end
    
    %                 rmdir(reduced_dir, 's')
    
                catch ME
                    fprintf(1,'The identifier was:\n%s', ME.identifier);
                    fprintf(1,'There was an error! The message was:\n%s', ME.message);
      
                    fid = fopen(gA.error_path,'a+');
                    fprintf(fid,'%s\n',ME.message);
                    for e=1:length(ME.stack)
                      fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                    end
                    
                    % close file
                    fclose(fid);
                end
            
            else

                for t_idx = time_points
                    
                    reduced_dir = [gA.segm_dir filesep 'reduced_images' filesep group_name filesep num2str(t_idx)];
                    gA.extractMaskCellpose(reduced_dir, raw_mask_dir)
        
                    try

        
        %                 rmdir(reduced_dir, 's')
        
                    catch ME
                        fprintf(1,'The identifier was:\n%s', ME.identifier);
                        fprintf(1,'There was an error! The message was:\n%s', ME.message);
          
                        fid = fopen(gA.error_path,'a+');
                        fprintf(fid,'%s\n',ME.message);
                        for e=1:length(ME.stack)
                          fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                        end
                        
                        % close file
                        fclose(fid);
                    end
                end

                parfor s_idx_idx = 1:length(indices)
                    s_idx = indices(s_idx_idx)

                    S1 = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                    final_dim = size(S1);
                    
                    path2label = [raw_mask_dir filesep gA.sample_names{s_idx} '*_cp_masks.tif'];
                    linkLabels(path2label, gA.segm_dir, [gA.sample_names{s_idx} '_mask.tif'], [], final_dim, gA.param.link_labels);
%                         delete(path2label)
                end
            end

        end
        
        %% analysis functions
        % Runs Morgana analysis to extract morphological and fluorescence properties
        function extractPropertiesMorgana(gA)

            gA.properties_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_properties'];
            mkdir(gA.properties_dir);
            
            try
                run_morgana_analysis_flourescence(gA.channel_dir, gA.segm_dir, gA.properties_dir, gA.param.plate_mode)
%                 run_morgana_analysis_flourescence(gA.channel_dir, gA.cellprob_segm_dir, gA.properties_dir, gA.param.plate_mode)
            
            catch ME
                fprintf(1,'The identifier was:\n%s', ME.identifier);
                fprintf(1,'There was an error! The message was:\n%s', ME.message);
  
                fid = fopen(gA.error_path,'a+');
                fprintf(fid,'%s\n',ME.message);
                for e=1:length(ME.stack)
                  fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                end
                
                % close file
                fclose(fid);
            end

        end
        
        % Aggregates analysis results from Morgana Json files by experimental groups, extracting shape, fluorescence, and spatial metrics
        function groupData(gA) % add coordinates of spine as cell arrays
            gA.group_analysis_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_group_analysis'];
            mkdir(gA.group_analysis_dir);

            try
                data_dir = gA.properties_dir;
                save_folder = gA.group_analysis_dir;
                        
                morph_json_file = [data_dir filesep 'aggregated_results' '_morpho_params.json'];
                fluo_json_file = [data_dir filesep 'aggregated_results' '_fluo_intensity.json'];
    
                fid = fopen(morph_json_file); 
                raw = fread(fid,inf); 
                str = char(raw'); 
                fclose(fid); 
                morg_res = jsondecode(str);
                
                fid = fopen(fluo_json_file); 
                raw = fread(fid,inf); 
                str = char(raw'); 
                fclose(fid); 
                fluo_res = jsondecode(str);
                
                %%
                
                results = struct();
                
                disp('Extracting data morgana...')
                % fill struct with morganas data
                counter = 1;
                dictionary = {};
                for page_id = 1:length(morg_res)
                    disp(page_id)
                    
                    % extract info from the original name file
                    getnum = @(x) str2double(x{1});
                    s_idx = getnum(extractBetween(morg_res(page_id).mask_file, filesep, '_'));

                    % infere experimental group
                    if gA.param.plate_mode
                        regex_res = cellfun(@(x) regexp(morg_res(page_id).mask_file,  ['[A-Z]' x '(?!\d)'])>0, gA.groups, 'UniformOutput', false); % determine group id
                        group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last');

                        get_first = @(x) x{1};
                        s_idx_idx = find(cell2mat(cellfun(@(x) str2double(get_first(strsplit(x,'_'))), gA.sample_names, 'UniformOutput', false))==s_idx);

                    
                    else
                        regex_res = cellfun(@(x) regexp(morg_res(page_id).mask_file, ['\d+_' x '_\d+'] )>0, gA.groups, 'UniformOutput', false); % determine group id
                        s_idx_idx = s_idx;
                        % regex_res = cellfun(@(x) extract(morg_res(page_id).mask_file, x ), gA.groups, 'UniformOutput', false); % determine group id
                        
                        group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last');
                        if isempty(group_id)
                            group_id = find(cellfun(@(x) contains(morg_res(page_id).mask_file, x), gA.groups));
                        end
                        % to test
                        % matches = regexp(gA.sample_names, ['^' num2str(s_idx)], 'once');
                        % s_idx_idx = ~cellfun('isempty', matches);
                        % s_idx_idx = find(startsWith(gA.sample_names,  {num2str(s_idx, '%02d'), num2str(s_idx, '%03d')}, 'IgnoreCase', true));

                    end

                    results(s_idx).group_id = group_id;
                    results(s_idx).group = gA.groups{group_id};
                    results(s_idx).name = gA.sample_names{s_idx_idx};
                    counter = counter+1;


                    % This block of code is because "A01" instead of "01_A01"
                    if isnan(s_idx)
                        wellname = extractBetween(morg_res(page_id).mask_file, '\', '.tif');
                        s_idx = find(cell2mat(cellfun(@(x) strcmp(wellname, x), dictionary, 'UniformOutput', false)));
                        if isempty(s_idx)
                            s_idx = counter;
                            dictionary{s_idx} = wellname{1};
                            
                            if gA.param.plate_mode
                                regex_res = cellfun(@(x) regexp(morg_res(page_id).mask_file,  ['[A-Z]' x '(?!\d)'])>0, gA.groups, 'UniformOutput', false); % determine group id
                                group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last'); 
                            end

                            results(s_idx).group_id = group_id;
                            results(s_idx).group = gA.groups{group_id};
                            results(s_idx).name = wellname{1};
                            counter = counter+1;
                        end
                    end


                    t_idx = getnum(extractBetween(morg_res(page_id).mask_file, '.tif_', '_cropped'));
                    g_idx = getnum(extractBetween(morg_res(page_id).mask_file, '_cropped', '_'));
                    
                    try

                        results(s_idx).sample(g_idx).eccentricity(t_idx) = morg_res(page_id).eccentricity(1); % E
                        results(s_idx).sample(g_idx).perimeter(t_idx) = morg_res(page_id).perimeter(1); % P
                        results(s_idx).sample(g_idx).major_axis_length(t_idx) = morg_res(page_id).major_axis_length(1); % MA
                        results(s_idx).sample(g_idx).minor_axis_length(t_idx) = morg_res(page_id).minor_axis_length(1); % ma
                        results(s_idx).sample(g_idx).axis_ratio(t_idx) = morg_res(page_id).major_axis_length(1)./morg_res(page_id).minor_axis_length(1); % ellipse ratio
                        results(s_idx).sample(g_idx).area(t_idx) = morg_res(page_id).area(1); % A
                        results(s_idx).sample(g_idx).form_factor(t_idx) = morg_res(page_id).form_factor(1); % FF
                        results(s_idx).sample(g_idx).shape_index(t_idx) = morg_res(page_id).perimeter(1)./sqrt(morg_res(page_id).area(1)); % shape index
                        results(s_idx).sample(g_idx).locoefa(:,t_idx) = morg_res(page_id).locoefa_coeff;

                        r = corrcoef( morg_res(page_id).midline(:,1),  morg_res(page_id).midline(:,2));
                        r2 = r(1,2)^2;

                        segments = diff(morg_res(page_id).midline);
                        path_length = sum(sqrt(sum(segments.^2, 2)));
                        direct_dist = sqrt(sum((morg_res(page_id).midline(end,:) - morg_res(page_id).midline(1,:)).^2));
                        path_efficiency = direct_dist / path_length;

                        results(s_idx).sample(g_idx).straightness_corr(t_idx) = r2;
                        results(s_idx).sample(g_idx).midline_length(t_idx) = path_length; % ML
                        results(s_idx).sample(g_idx).path_efficiency(t_idx) = path_efficiency;
                        results(s_idx).sample(g_idx).lengthening_ratio(t_idx) = path_length./results(s_idx).sample(g_idx).midline_length(1); % lengthening ratio
                    catch
                        disp(['Properties not computed for: ' morg_res(page_id).mask_file])

                        results(s_idx).sample(g_idx).eccentricity(t_idx) = nan; % E
                        results(s_idx).sample(g_idx).perimeter(t_idx) = nan; % P
                        results(s_idx).sample(g_idx).major_axis_length(t_idx) = nan; % MA
                        results(s_idx).sample(g_idx).minor_axis_length(t_idx) = nan; % ma
                        results(s_idx).sample(g_idx).axis_ratio(t_idx) = nan; % ellipse ratio
                        results(s_idx).sample(g_idx).area(t_idx) = nan;% A
                        results(s_idx).sample(g_idx).form_factor(t_idx) = nan; % FF
                        results(s_idx).sample(g_idx).shape_index(t_idx) = nan; % shape index
                        results(s_idx).sample(g_idx).locoefa(:,t_idx)  = nan;

                        results(s_idx).sample(g_idx).straightness_corr(t_idx) = nan;% ML
                        results(s_idx).sample(g_idx).midline_length(t_idx) = nan;
                        results(s_idx).sample(g_idx).path_efficiency(t_idx) = nan;
                        results(s_idx).sample(g_idx).lengthening_ratio(t_idx)  = nan; % lengthening ratio
                    end
                    
                    
                    % security check
                    s_idx = getnum(extractBetween(fluo_res(page_id).mask_file, filesep, '_'));
                    t_idx = getnum(extractBetween(fluo_res(page_id).mask_file, '.tif_', '_cropped'));
                    g_idx = getnum(extractBetween(fluo_res(page_id).mask_file, '_cropped', '_'));
                
                    % extract fluorescence data
                    for ch_idx = 0:gA.param.n_channels-1
                        try
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_ANGprofile']){t_idx} = fluo_res(page_id).(['ch' num2str(ch_idx) '_ANGprofile']);
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_APprofile']){t_idx}  = fluo_res(page_id).(['ch' num2str(ch_idx) '_APprofile']);
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_Average']){t_idx}    = fluo_res(page_id).(['ch' num2str(ch_idx) '_Average']);
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_Background']){t_idx} = fluo_res(page_id).(['ch' num2str(ch_idx) '_Background']);
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_LRprofile']){t_idx}  = fluo_res(page_id).(['ch' num2str(ch_idx) '_LRprofile']);
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_RADprofile']){t_idx} = fluo_res(page_id).(['ch' num2str(ch_idx) '_RADprofile']);
    
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx)])(t_idx) = [fluo_res(page_id).(['ch' num2str(ch_idx) '_Average']) - fluo_res(page_id).(['ch' num2str(ch_idx) '_Background'])];
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_AP_pol'])(t_idx) = mean(sqrt(diff(movmean([fluo_res(page_id).(['ch' num2str(ch_idx) '_APprofile']) - fluo_res(page_id).(['ch' num2str(ch_idx) '_Background'])],100)).^2));
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_LR_pol'])(t_idx) = mean(sqrt(diff(movmean([fluo_res(page_id).(['ch' num2str(ch_idx) '_LRprofile']) - fluo_res(page_id).(['ch' num2str(ch_idx) '_Background'])],100)).^2));
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_RD_pol'])(t_idx) = mean(sqrt(diff(movmean([fluo_res(page_id).(['ch' num2str(ch_idx) '_RADprofile']) - fluo_res(page_id).(['ch' num2str(ch_idx) '_Background'])],100)).^2));
                        
                        catch
%                             disp(['Properties not computed for: ' morg_res(page_id).mask_file])

                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_ANGprofile']){t_idx} = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_APprofile']){t_idx}  = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_Average']){t_idx}    = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_Background']){t_idx} = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_LRprofile']){t_idx}  = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_RADprofile']){t_idx} = nan;
    
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx)])(t_idx) = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_AP_pol'])(t_idx) = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_LR_pol'])(t_idx) = nan;
                            results(s_idx).sample(g_idx).(['ch' num2str(ch_idx) '_RD_pol'])(t_idx) = nan;
                        end
                    
                    end

                end
                
                %%
                if ~isempty(gA.semant_dir) % deprecated
                    disp('Extracting data semantic segmentation...')
                    for s_idx_idx = 1:length(gA.sample_names) % par for?
                        s_idx  = str2double(extractBefore(gA.sample_names{s_idx_idx},'_'));
                        % security check
                        if isnan(s_idx)
                            wellname = gA.sample_names{s_idx_idx};
                            s_idx = find(cell2mat(cellfun(@(x) strcmp(wellname, x), dictionary, 'UniformOutput', false)));
                        end
%                             
                        % open stacks
                        % S = openTiffStack([stack_list(s_idx).folder filesep stack_list(s_idx).name]);
                        O = openTiffStack([gA.segm_dir filesep gA.sample_names{s_idx_idx} '_mask.tif']);
                        B = openTiffStack([gA.semant_dir filesep gA.sample_names{s_idx_idx} '_mask.tif']);
                        
                        % determine number of labels within data set
                        g_inds = unique(O(O>0));
                    
                        % extract the data
                        for t_idx = 1:size(O, 3)
                            
                            Ci = B(:,:,t_idx)==255;
                            for g_idx = g_inds'
                    
                                     % extract spreading area
                                    Mi = imdilate(O(:,:,t_idx)==g_idx, strel('disk', 5));
                                    props = regionprops(Mi, 'Centroid', 'Circularity');
                                    
                                    if length(props) == 1
                        
                                        % cell and body data
                                        results(s_idx).sample(g_idx).coords(t_idx, :) = props.Centroid; % XY
                                        results(s_idx).sample(g_idx).area_body(t_idx) = sum(Mi(:)); % area total
                                        results(s_idx).sample(g_idx).area_spreading(t_idx)  = sum(Mi(:) & Ci(:)); % area cells
                                        results(s_idx).sample(g_idx).circularity(t_idx)  = props.Circularity; % circularity
                                        results(s_idx).sample(g_idx).spreading_ratio(t_idx) = sum(Mi(:) & Ci(:))./sum(Mi(:)); % spreading
                                    end
                    
                            end
                        end
                    
                    end
                end

                if ~isempty(gA.cellprob_segm_dir)    
                    disp('Extracting data cell probability segmentation...')
                    for s_idx_idx = 1:length(gA.sample_names) % par for?
                        
                        s_idx  = str2double(extractBefore(gA.sample_names{s_idx_idx},'_'));
                        
                        % security check
                        if isnan(s_idx)
                            wellname = gA.sample_names{s_idx_idx};
                            s_idx = find(cell2mat(cellfun(@(x) strcmp(wellname, x), dictionary, 'UniformOutput', false)));
                        end
                            
                        % open stacks
                        MB = openTiffStack([gA.segm_dir filesep gA.sample_names{s_idx_idx} '_mask.tif']); % mask body
%                         MC = openTiffStack([gA.cellprob_segm_dir filesep gA.sample_names{s_idx_idx} '_mask.tif']); % mask cells
                       

                        % determine number of labels within data set
                        g_inds = unique(MB(MB>0));
                    

                        % extract the data
                        for t_idx = 1:size(MB, 3)
                            
                            MC_t_idx = zeros(size(MB, 1), size(MB, 2));
                            try 
                                mask_filename = [gA.cellprob_segm_dir filesep 'raw_mask' filesep gA.sample_names{s_idx_idx} '_t' num2str(t_idx, '%04d') '.tif_mask.tif'];
                                MC_t_idx = imread(mask_filename);
                                MC_t_idx = imresize(MC_t_idx, [size(MB, 1) size(MB, 2)], 'nearest');
                            catch
                                disp([gA.sample_names{s_idx_idx} '_t' num2str(t_idx, '%04d') '.tif_mask.tif' ' not computed'])
                            end
                            Ci_all = double(MC_t_idx==1);

                            % dist transform to split cell mask
                            dist_t = zeros(size(MB,1), size(MB,2), length(g_inds));
                            for g_idx = g_inds'
                                dist_t(:,:,g_idx) = bwdist(MB(:,:,t_idx)==g_idx);
                            end
                            [~, ind_map] = min(dist_t, [], 3);
                            
                            for g_idx = g_inds'
                                    % extact the part of the image corresponding to an index
                                    Ci = Ci_all.*ind_map==g_idx;
                                    ci_props = regionprops(Ci, 'Circularity');
                                    if isempty(ci_props)
                                        ci_props(1).Circularity = nan;
                                    end

                                    % extract spreading area
                                    Mi = double(MB(:,:,t_idx)==g_idx);
                                    props = regionprops(Mi, 'Centroid', 'Circularity');
                                    
                                    if length(props) == 1
                        
                                        % cell and body data
                                        results(s_idx).sample(g_idx).coords(t_idx, :) = props.Centroid; % XY
                                        results(s_idx).sample(g_idx).area_body(t_idx) = sum(Mi(:)); % area total
                                        results(s_idx).sample(g_idx).area_spreading(t_idx)  = sum(Ci(:) & ~Mi(:)); % area cells
                                        results(s_idx).sample(g_idx).circularity(t_idx)  = props.Circularity; % circularity
                                        results(s_idx).sample(g_idx).spreading_ratio(t_idx) = sum(Ci(:) & ~Mi(:))./( sum(Ci(:) & ~Mi(:)) + sum(Mi(:))); % spreading
                                        results(s_idx).sample(g_idx).circularity_spreading(t_idx) = ci_props.Circularity; % shape_index_spreading
                                    end
                    
                            end
                        end
                    
                    end
                end

                
                save([save_folder filesep 'results'], 'results');

            catch ME
                fprintf(1,'The identifier was:\n%s', ME.identifier);
                fprintf(1,'There was an error! The message was:\n%s', ME.message);
  
                fid = fopen(gA.error_path,'a+');
                fprintf(fid,'%s\n',ME.message);
                for e=1:length(ME.stack)
                  fprintf(fid,'%sin %s at %i\n',ME.message,ME.stack(e).name,ME.stack(e).line);
                end
                
                % close file
                fclose(fid);
            end

        end

        function resegmentData(gA, plate_indices, new_diameter) % TO FINISH
            present_plate_indices = cellfun(@(x) sscanf(x,'%d'), gA.sample_names);
            experiment_indices = find(ismember(present_plate_indices, plate_indices));
            present_plate_indices(experiment_indices) % to check
            
            gA.current_idx = 3;
%             gA.findGastruloids([], new_diameter, experiment_indices);
            gA.findGastruloidBodies([], new_diameter, experiment_indices);

            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);

            parfor t_idx = 1:size(S, 3)
                if gA.param.plate_mode
                    gA.plotSummaryMasks96WellPlate(t_idx, 0.1);
                else
                    gA.plotSummaryMasks(t_idx, 0.1);
                end
            end
            
            gA.current_idx = 5;
            gA.extractPropertiesMorgana();
            gA.stepIdx(); gA.save()
            
            %%
            gA.groupData();
            gA.save()
            
            gA.plotSummary()
        end

        
        %% Visualization Methods
        % Creates box plots comparing properties across groups at specified timepoint
        function plotPropertiesBoxPlots(gA, t_idx, plot_prop_names, remove_outliers)
            %% boxplot final state
            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable
            
            param = gA.group_analysis_param;
            save_folder = gA.group_analysis_dir;
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            
            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end
            
            disp('Plotting bloxplots...')
            
            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = size(S, 3); end
            center = [size(S,1), size(S,2)]./2;
            
            
                      
           if ~exist('plot_prop_names', 'var') || isempty(plot_prop_names)
                plot_prop_names = {'area', 'spreading_ratio', 'circularity',...
                                   'ch1', 'ch2','ch1_AP_pol', 'ch1_RD_pol', 'ch2_RD_pol'};
                
            end        
            
            if ~exist('remove_outliers', 'var') || isempty(remove_outliers);
                remove_outliers = 1;
            end


            for prop = plot_prop_names
               
               group_data = nan(length(gA.groups), 100000);
               group_counter = ones(1,length(gA.groups));
               for s_idx = 1:length(results)
                    group_id = results(s_idx).group_id;
            
                   for g_idx = 1:length(results(s_idx).sample)
                            % this is for ~plate mode
                            %C = data(ri).data(gi,2:3,1); % centroid at t=0
                            %A = data(ri).data(gi,4,1); % Area at t=0
                            %dist = sqrt((C(1) - center(1)).^2 + (C(2) - center(2)).^2);
                            try
                                time_data = results(s_idx).sample(g_idx).(prop{1});
                                if length(time_data) >= t_idx
                                    if ~isnan(time_data(t_idx))
                                        group_data(group_id, group_counter(group_id)) = time_data(t_idx);
                                        group_counter(group_id) = group_counter(group_id)+1;
                                    end
                                end
                            end   
                            
                   end
               end
            
                group_data(:,max(group_counter):end) = [];

                clf
                boxplot(group_data', gA.groups,...
                        'ColorGroup', [], 'Colors', gA.param.cmap,...
                        'BoxStyle', 'outline', 'Symbol', 'k+')
                xticklabels(gA.group_final_names)
                
                hold on
                for gi = 1:length(gA.groups)
                    data2plot = group_data(gi,:);

                    plot(gi+rand(1,length(data2plot))./5-0.5/5, data2plot, ...
                        '.', 'MarkerSize', 20, 'Color', cmap(gi,:));
                    
                    xlim([0 length(gA.groups)+1])
                end
                

                if remove_outliers
                    [lim_data, ~]= rmoutliers(group_data', 'gesd');
                end

                try
                    ylim([min(lim_data(:)), max(lim_data(:))])
                end
                
                for ch_idx = 0:gA.param.n_channels-1
                    prop = strrep(prop, ['ch' num2str(ch_idx)], gA.param.channel_names{ch_idx+1});
                end
            
                box on; grid on;
                ylabel(strrep(prop, '_', ' '))
                
                set(gcf, 'Color', 'w')
                fontsize(gcf, 14, 'points')
            
                export_fig([save_folder filesep 'boxplot_final_' prop{1} '_t' num2str(t_idx, '%02d') '_.png'])
            end

        end

        % Generates PCA plots of morphological/fluorescence properties
        function plotPropertiesPCA(gA, t_plot, remove_outliers, plot_prop_names)
            % use t_plot = -1 to PCA all variables at all time points
            %% PCA final state
            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable
            

            param = gA.group_analysis_param;
            save_folder = gA.group_analysis_dir;
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            
            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end
            
            if ~exist('t_plot', 'var') || isempty(t_plot)
                t_plot = size(S,3);
            end

            if ~exist('remove_outliers', 'var') || isempty(remove_outliers)
                remove_outliers = false;
            end

            
            center = [size(S,1), size(S,2)]./2;

            group_data = nan(length(gA.groups), 100000);
            
            if ~exist('plot_prop_names', 'var') || isempty(plot_prop_names)
                plot_prop_names = {'area', 'spreading_ratio', 'circularity',...
                                   'ch1', 'shape_index'};
                
            end            
            
            data_matrix = zeros(0,length(plot_prop_names)*size(S,3));
            property_id = 1;
            group_indices = zeros(0,1);
            sample_indices = zeros(0,1);
            gastruloid_indices = zeros(0,1);
            for prop = plot_prop_names
               
               sample_counter = 1;
               for s_idx = 1:length(results)
                    group_id = results(s_idx).group_id;
                    group_id
                   for g_idx = 1:length(results(s_idx).sample)
                            % this is for ~plate mode
                            %C = data(ri).data(gi,2:3,1); % centroid at t=0
                            %A = data(ri).data(gi,4,1); % Area at t=0
                            %dist = sqrt((C(1) - center(1)).^2 + (C(2) - center(2)).^2);
                            try
                                data_matrix(sample_counter, property_id+[1:length(results(s_idx).sample(g_idx).area)]-1) = nan;
                                for t_idx = 1:length(results(s_idx).sample(g_idx).area)
                                    time_data = results(s_idx).sample(g_idx).(prop{1});
                                    
                                    if ~isnan(time_data(t_idx))
                                        data_matrix(sample_counter, property_id+t_idx-1) = time_data(t_idx);
                                    end
    
                                    
                                end
                                
                                group_indices(sample_counter) = group_id;
                                sample_indices(sample_counter) = s_idx;
                                gastruloid_indices(sample_counter) = g_idx;
                                sample_counter = sample_counter + 1;
                            end
                   end
                   
               end
               property_id = property_id+size(S,3);
            end

            size(group_indices)

            size(data_matrix)
            %
            disp('Plotting PCA...')
            X = (data_matrix-mean(data_matrix, "omitnan"))./std(data_matrix, "omitnan");
            X(isnan(X(:))) = 0;
%             kk = [11, 60, 65, 86];
%             X(kk,:) = [];
%             group_indices(kk) = [];

            
            if remove_outliers
                final_group_indices = [];
                final_X = zeros(0, size(X,2));
                for g_id = 1:length(gA.groups)
                    g_indices = find(group_indices==g_id);
                    [gX, rm_inds]= rmoutliers(X(g_indices,:), 'gesd', 'ThresholdFactor', 0.5); %very aggresive 0.5
                    g_indices(rm_inds) = [];

                    final_X(end+1:end+size(gX, 1),:) = gX;
                    final_group_indices(end+1:end+size(gX, 1)) = group_indices(g_indices);
                end
                X = final_X;
                group_indices = final_group_indices;

            end
            
            if t_plot ~= -1
                X = X(:,t_plot:size(S,3):end);
                plot_prop_names_time = plot_prop_names;
            else
                plot_prop_names_time = {};
                for t_idx = 1:size(S,3)
                    for pn = plot_prop_names
                        plot_prop_names_time{end+1} = [pn{1} '_' num2str(t_idx)];
                    end        
                end
            end

            for prop_id = 1:length(plot_prop_names_time)
                for ch_idx = 0:gA.param.n_channels-1
                    plot_prop_names_time{prop_id} = strrep(plot_prop_names_time{prop_id}, ['ch' num2str(ch_idx)], gA.param.channel_names{ch_idx+1});
                end
            end

            clf
%             [coeff,score,latent] = ppca(X,2);
            [coeff,score,latent] = pca(X);

%             if remove_outliers
%                 [score, rm_inds]= rmoutliers(score, 'percentile', [20, 80]);
%                 group_indices(rm_inds) = [];
%             end

            if t_plot ~= -1
                biplot(coeff(:,1:2), 'scores', score(:,1:2), 'VarLabels', plot_prop_names_time, 'LineWidth', 2, 'Color', 'k');
            else
                 biplot(coeff(:,1:2), 'scores', score(:,1:2), 'LineWidth', 2, 'Color', 'k');
            end
            
            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end
            
            h = findobj(gca,"Tag","obsmarker");
            

            hold on
            scatter(h.XData, h.YData, 60, cmap(group_indices,:), 'filled', 'o')

            for g_id = 1:length(gA.groups)
                hold on
                plot(median(h(1).XData(group_indices==g_id)), median(h(1).YData(group_indices==g_id)), 'sq', 'MarkerFaceColor', cmap(g_id,:), 'MarkerEdgeColor','k', 'MarkerSize', 12)
                points = [h(1).XData(group_indices==g_id); h(1).YData(group_indices==g_id)]';
                size(points)
                scale_factor = 1;
%                 [center_pca, axes_pca, rot_pca] = fitEnclosingEllipse(points, 'pca', scale_factor);
%                 plotEllipse(points, center_pca, axes_pca, rot_pca, 'pca');

                [center_pca, axes_pca, rot_pca] = fitEnclosingEllipse(points, 'mvee', scale_factor);
                he = plotEllipse(points, center_pca, axes_pca, rot_pca, 'mvee');
                he.Color = cmap(g_id,:);
                he.LineWidth = 1;
            end
            
            legend_cmap = cmap;
            qw = {};
            legend_names = gA.group_final_names;
            for idx = 1:length(gA.groups)
                legend_names{idx} = strjoin(strsplit(legend_names{idx},'_'), ' ');
                qw{idx} = scatter(nan, nan, 40, legend_cmap(idx,:), 'filled', 'o');
            end

            legend([qw{:}], legend_names, 'location', 'bestoutside', 'box', 'off');
            grid on; box on;
            set(gcf, 'Color', 'w')
            fontsize(gcf, 14, 'points')
            title(num2str(t_plot, '%02d'))


                        
            export_fig([save_folder filesep 'PCA_' num2str(t_plot, '%02d') '.png']);

        end
        
        % Plots PCA trajectories showing temporal evolution of properties
        function plotPCATrajectories(gA, remove_outliers, plot_prop_names)
            %% PCA trajectories
            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable
            
            param = gA.group_analysis_param;
            save_folder = gA.group_analysis_dir;
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            
            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end
            
            if ~exist('remove_outliers', 'var') || isempty(remove_outliers)
                remove_outliers = false;
            end
            
            center = [size(S,1), size(S,2)]./2;

            group_data = nan(length(gA.groups), 100000);
            

%             plot_prop_names = {'eccentricity', 'perimeter', 'major_axis_length', ...
%             'minor_axis_length', 'axis_ratio', 'perimeter', 'midline_length',...
%             'area', 'form_factor', 'shape_index', 'lengthening_ratio',...
%             'bra', 'sox', 'bra_pol', 'sox_pol',...
%             'area_body', 'area_spreading', 'circularity', 'spreading_ratio'};

            if ~exist('plot_prop_names', 'var') || isempty(plot_prop_names)
                plot_prop_names = {'area', 'spreading_ratio', 'circularity',...
                                   'ch1', 'ch2','ch1_AP_pol', 'ch1_RD_pol', 'ch2_RD_pol'}; 
            end       
            
            data_matrix = zeros(0,length(plot_prop_names));
            property_id = 1;
            time_indices = zeros(0,1);
            group_indices = zeros(0,1);
            for prop = plot_prop_names
               
               sample_counter = 1;
               for s_idx = 1:length(results)
                    group_id = results(s_idx).group_id;
            
                    
                   for g_idx = 1:length(results(s_idx).sample)
                            % this is for ~plate mode
%                             C = results(s_idx).sample(g_idx).('coords');
%                             C = C(1,:); % centroid at t=0
%                             
%                             A = results(s_idx).sample(g_idx).('area'); 
%                             A = % Area at t=0
                            %dist = sqrt((C(1) - center(1)).^2 + (C(2) - center(2)).^2);
                            for t_idx = 1:length(results(s_idx).sample(g_idx).area)
                                time_data = results(s_idx).sample(g_idx).(prop{1});
                                if ~isnan(time_data(t_idx)) && time_data(t_idx)~=0
                                    data_matrix(sample_counter, property_id) = time_data(t_idx);
                                    time_indices(sample_counter) = t_idx;
                                    group_indices(sample_counter) = group_id;
                                    sample_counter = sample_counter + 1;
                                end

                                
                            end
                   end
                   
               end
               property_id = property_id+1;
            end
            
            %%
            X = (data_matrix-nanmean(data_matrix))./nanstd(data_matrix);
            X(isnan(X(:))) = 0;
%             [X, inds] = unique(X, 'rows');
%             time_indices = time_indices(inds);
%             group_indices = group_indices(inds);

            if remove_outliers
                [X, rm_inds]= rmoutliers(X, 'median');
                group_indices(rm_inds) = [];
                time_indices(rm_inds) = [];
            end

            clf
            [coeff,score,latent] = pca(X);
            h = {};
            
            figure(1); clf;

            trajs = zeros(length(gA.groups), size(S,3), 2);
            for g_ind = 1:length(gA.groups)
                for t_idx = 1:size(S,3)
                    current_ids = (time_indices==t_idx) & (group_indices==g_ind);
                    hold on
                    figure(2)
                    if g_ind ==1 && t_idx ==1; figure(1); end

                    if remove_outliers
                        [scores, rm_inds]= rmoutliers(score(current_ids,1:2), 'median');
                        current_ids(rm_inds) = [];
                    else
                        scores = score(current_ids,1:2);
                    end

                    for prop_id = 1:length(plot_prop_names)
                        for ch_idx = 0:gA.param.n_channels-1
                            plot_prop_names{prop_id} = strrep(plot_prop_names{prop_id}, ['ch' num2str(ch_idx)], gA.param.channel_names{ch_idx+1});
                        end
                    end

                    biplot(coeff(:,1:2), 'scores', score(current_ids,1:2), 'VarLabels', plot_prop_names, 'LineWidth', 0.1, 'Color', 'k', 'Marker','none');
                    
                    h = findobj(gca,"Tag","obsmarker");
        
                    figure(1)
                    hold on
                    scatter(h(1).XData, h(1).YData, 60, cmap(group_indices(current_ids),:), 'filled', 'o', 'MarkerFaceAlpha', 0.2)
    
                    trajs(g_ind, t_idx, 1) = median(h(1).XData);
                    trajs(g_ind, t_idx, 2) = median(h(1).YData);

                    grid on; box on;
                    set(gcf, 'Color', 'w')
                    fontsize(gcf, 14, 'points')
                end
            end

            
            for g_ind = 1:length(gA.groups)
                hold on
                plot(trajs(g_ind,1,1), trajs(g_ind,1,2), 'o', 'color', cmap(g_ind,:), 'MarkerSize', 12, 'LineWidth', 2)
                plot(trajs(g_ind,end,1), trajs(g_ind,end,2), 'x', 'color', cmap(g_ind,:), 'MarkerSize', 12, 'LineWidth', 2)
                plot(trajs(g_ind,:,1), trajs(g_ind,:,2), '.-', 'color', cmap(g_ind,:), 'MarkerSize', 12, 'LineWidth', 2)
            end

            legend_cmap = cmap;
            qw = {};
            legend_names = gA.group_final_names;
            for idx = 1:length(gA.groups)
                legend_names{idx} = strjoin(strsplit(legend_names{idx},'_'), ' ');
                qw{idx} = scatter(nan, nan, 40, legend_cmap(idx,:), 'filled', 'o');
            end

            legend([qw{:}], legend_names, 'location', 'bestoutside');

            export_fig([save_folder filesep 'PCAtrajectories.png']);            


        end

        % Generates comprehensive analysis plots including PCA, box plots, and time series
        function plotSummary(gA, plot_prop_names)

            % plot all PCA plots removing outliers
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            remove_outliers = 0;

            if ~exist('plot_prop_names', 'var') || isempty(plot_prop_names)
                plot_prop_names = {'area', 'spreading_ratio', 'shape_index', 'circularity', 'midline_length',...
                           'ch1', 'ch2', 'area_spreading', 'straightness_corr', 'path_efficiency'};
            end      
                % plot_prop_names = {'area_body', 'spreading_ratio', 'shape_index', 'circularity', 'midline_length',...
                %        'ch1', 'ch2','ch1_AP_pol', 'ch1_RD_pol', 'ch2_AP_pol', 'ch2_RD_pol', 'area_spreading'};

%             plot_prop_names = {'area', 'shape_index', 'midline_length', 'ch1', 'axis_ratio'};
            for ch_idx = 1:gA.param.n_channels-1 % exclude first channel, always brigthfield
                plot_prop_names{end+1} = ['ch' num2str(ch_idx)];
                plot_prop_names{end+1} = ['ch' num2str(ch_idx) '_AP_pol'];
            end

            
            for t_idx = 1:size(S, 3)
                gA.plotPropertiesPCA(t_idx, remove_outliers, plot_prop_names);
                gA.plotPropertiesBoxPlots(t_idx, plot_prop_names);
            end
            
            gA.plotPropertiesPCA(-1, remove_outliers, plot_prop_names); %all together
%             gA.plotPCATrajectories(1, plot_prop_names);
            if size(S,3) > 1; gA.plotTimeLines(plot_prop_names); end
              
        end
        
        %  Creates time-course plots showing property evolution over time
        function plotTimeLines(gA, plot_prop_names, remove_outliers)
            
            %% line plots
            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable
            
            param = gA.group_analysis_param;
            save_folder = gA.group_analysis_dir;
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            
            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end
            
            if ~exist('t_plot', 'var') || isempty(t_plot)
                t_plot = size(S,3);
            end

            if ~exist('remove_outliers', 'var') || isempty(remove_outliers)
                remove_outliers = false;
            end

            
            center = [size(S,1), size(S,2)]./2;

            
            
            if ~exist('plot_prop_names', 'var') || isempty(plot_prop_names)
                plot_prop_names = {'area', 'spreading_ratio', 'circularity',...
                                   'ch1', 'ch2','ch1_AP_pol', 'ch1_RD_pol', 'ch2_RD_pol'};
            end            
            
            [rr, cc] = number_subplots(length(gA.groups));      
            
            
            for prop = plot_prop_names
               clf
               group_indices = zeros(0,1);
               group_data = nan(length(gA.groups), 100000, size(S,3));
               sample_counter = 1;
               for s_idx = 1:length(results)
                    group_id = results(s_idx).group_id;
            
                   for g_idx = 1:length(results(s_idx).sample)
                            % this is for ~plate mode
                            %C = data(ri).data(gi,2:3,1); % centroid at t=0
                            %A = data(ri).data(gi,4,1); % Area at t=0
                            %dist = sqrt((C(1) - center(1)).^2 + (C(2) - center(2)).^2);
                            for t_idx = 1:length(results(s_idx).sample(g_idx).area)
                                time_data = results(s_idx).sample(g_idx).(prop{1});
                                
                                subplot(rr, cc, group_id)
                                time_data = movmean(time_data, 1);

                                hold on
                                plot(time_data, ...
                                    'o-',...
                                    'color', cmap(group_id,:).*(rand()*0.4+0.4), ...
                                    'MarkerSize', 2, 'LineWidth',0.1);
            
                                group_data(group_id, sample_counter, 1:length(time_data)) = time_data;

                                
                            end
                            group_indices(sample_counter) = group_id;
                            sample_counter = sample_counter + 1;
                   end
                   
               end

               group_data(group_data==0) = nan;
               group_data(:,sample_counter:end,:) = [];
                
                for gi = 1:length(gA.groups)
                    subplot(rr, cc, gi)
                    hold on
                    t_data = nanmedian(squeeze(group_data(gi,:,:)));
                    plot(t_data,...
                        'color', cmap(gi,:), 'LineWidth', 4)
                    if remove_outliers
                        [lim_data, ~]= rmoutliers(squeeze(group_data(gi,:,:)), 'movmedian', 3);
                        ylim([min(lim_data(:)), max(lim_data(:))])
                    end
                end
            
                % extract the limits
                for si = 1:length(gA.groups)
                    subplot(rr, cc, si)
                    if si == 1
                        xl = xlim;
                        yl = ylim;
                    elseif si > 1
                        xls = xlim;
                        yls = ylim;
                        if xls(1)<xl(1); xl(1) = xls(1); end
                        if xls(2)>xl(2); xl(2) = xls(2); end
                        if yls(1)<yl(1); yl(1) = yls(1); end
                        if yls(2)>yl(2); yl(2) = yls(2); end
                    end
                end
            
                for ch_idx = 0:gA.param.n_channels-1
                    prop = strrep(prop, ['ch' num2str(ch_idx)], gA.param.channel_names{ch_idx+1});
                end

                % set same axis
                for si = 1:length(gA.groups)
                    subplot(rr, cc, si)
                    if ~isempty(gA.group_final_names)
                        title(strjoin(strsplit(gA.group_final_names{si},'_'), ' '))
                    else
                        title(gA.groups{si},'_')
                    end

                    box on; grid on;
                    ylabel(strrep(prop, '_', ' '))
                    xlabel('frame')
                    xlim(xl)
                    ylim(yl)
                end
               
                set(gcf, 'Color', 'w')
                fontsize(gcf, 14, 'points')
                export_fig([save_folder filesep 'time_line_' prop{1} '.png'])
            end
        end

        % PCA analysis of LocoEFA shape coefficients with multiple normalization methods
        function plotLocoEfaPCA(gA, t_idx, first_mode, last_mode, group_ids)

            if ~exist('group_ids', 'var') || isempty(group_ids)
                group_ids = 1:length(gA.groups);
            end

            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable

            clf
            k = [];
            group_indices = [];
            area_norm = [];
            c = 1;
            for i = 1:length(results)
                if ~isempty(results(i).sample) & ismember(results(i).group_id, group_ids)
                    try
                        k(c,:) = results(i).sample.locoefa(:, t_idx);
                        group_indices(c) = results(i).group_id;
                        area_norm(c) = sqrt(results(i).sample.area(t_idx));
                        c = c+1;
                    end
                end
            end
            
            cmap = gA.param.cmap;
            if ~exist('first_mode', 'var') || isempty(first_mode)
                first_mode = 1;
            end


            if ~exist('last_mode', 'var') || isempty(last_mode)
                last_mode = 25;
            end
            
            t = tiledlayout(2, 2, 'TileSpacing','tight');
            title(t, ['Modes: ' num2str(first_mode) '-' num2str(last_mode)])
            
            nexttile
            data_matrix = k(:,first_mode+1:last_mode+1);
            X = data_matrix;
            X(isnan(X(:))) = 0;
            [coeff,score,latent] = pca(X);
            biplot(coeff(:,1:2), 'scores', score(:,1:2), 'LineWidth', 2, 'Color', 'k', 'VarLabels', num2str([first_mode:last_mode]'));
            h = findobj(gca,"Tag","obsmarker");
            
            hold on
            scatter(h.XData, h.YData, 60, cmap(group_indices,:), 'filled', 'o');
            title('RAW DATA')
            box on
            
            nexttile
            data_matrix = k(:,first_mode+1:last_mode+1);
            X = (data_matrix-mean(data_matrix, "omitnan"))./std(data_matrix, "omitnan");
            X(isnan(X(:))) = 0;
            [coeff,score,latent] = pca(X);
            biplot(coeff(:,1:2), 'scores', score(:,1:2), 'LineWidth', 2, 'Color', 'k', 'VarLabels', num2str([first_mode:last_mode]'));
            h = findobj(gca,"Tag","obsmarker");
            
            hold on
            scatter(h.XData, h.YData, 60, cmap(group_indices,:), 'filled', 'o');
            title('Z NORMALISED DATA')
            box on
            
            nexttile
            data_matrix = k(:,first_mode+1:last_mode+1);
            data_matrix = bsxfun(@rdivide, data_matrix, area_norm(:));
            X = data_matrix;
            X(isnan(X(:))) = 0;
            [coeff,score,latent] = pca(X);
            biplot(coeff(:,1:2), 'scores', score(:,1:2), 'LineWidth', 2, 'Color', 'k', 'VarLabels', num2str([first_mode:last_mode]'));
            h = findobj(gca,"Tag","obsmarker");
            
            hold on
            scatter(h.XData, h.YData, 60, cmap(group_indices,:), 'filled', 'o');
            title('AREA NORMALISED DATA')
            box on
            
            
            nexttile
            data_matrix = k(:,first_mode+1:last_mode+1);
            data_matrix = bsxfun(@rdivide, data_matrix, area_norm(:));
            X = (data_matrix-mean(data_matrix, "omitnan"))./std(data_matrix, "omitnan");
            X(isnan(X(:))) = 0;
            [coeff,score,latent] = pca(X);
            biplot(coeff(:,1:2), 'scores', score(:,1:2), 'LineWidth', 2, 'Color', 'k', 'VarLabels', num2str([first_mode:last_mode]'));
            h = findobj(gca,"Tag","obsmarker");
            
            hold on
            scatter(h.XData, h.YData, 60, cmap(group_indices,:), 'filled', 'o');
            title('AREA & Z NORMALISED DATA')
            box on

            legend_cmap = gA.param.cmap;
            qw = {};
            legend_names = gA.group_final_names;
            for idx = 1:length(gA.groups)
                legend_names{idx} = strjoin(strsplit(legend_names{idx},'_'), ' ');
                qw{idx} = scatter(nan, nan, 40, legend_cmap(idx,:), 'filled', 'o');
            end

            legend([qw{:}], legend_names, 'location', 'bestoutside');

%             set(h(:), {'Color'}, mat2cell(colors,ones(12,1),3))
            lg  = legend([qw{:}], legend_names,'Orientation','Horizontal','NumColumns',3); 
            lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

            set(gcf, 'Color', 'w')

            export_fig([gA.group_analysis_dir filesep 'LOCOEFA_PCA' num2str(t_idx, '%02d') '_' num2str(first_mode) '_' num2str(last_mode) '.png']);

        end
        
        % Analyzes shape complexity using higher-order LocoEFA modes
        function plotLocoEfaComplexity(gA, t_idx)
            % add modes from the 3rd onwards
            load([gA.group_analysis_dir filesep 'results.mat']); % loads results variable

            k = [];
            group_indices = [];
            area_norm = [];
            c = 1;
            group_data = nan(length(gA.groups), 100000);
            for i = 1:length(results)
                if ~isempty(results(i).sample)
                    try
                        k(c,:) = results(i).sample.locoefa(:, t_idx);
                        group_indices(c) = results(i).group_id;
                        area_norm(c) = sqrt(results(i).sample.area(t_idx));
                        c = c+1;
                    end
                end
            end
            
            data_matrix = k(:,2:51);
            data_matrix = bsxfun(@rdivide, data_matrix, area_norm(:));
            cmap = gA.param.cmap;
            
            first_mode_idx = 1;

            group_data = nan(length(gA.groups), 100000);
            group_counter = ones(length(gA.groups), 1);
            for s_idx = 1:length(group_indices)
                group_id =group_indices(s_idx);
                group_data(group_id, group_counter(group_id)) = sum(data_matrix(s_idx,first_mode_idx:end),2);
                group_counter(group_id) = group_counter(group_id)+1;
            end

            clf
            boxplot(group_data', gA.group_final_names,...
                    'ColorGroup', [], 'Colors', gA.param.cmap,...
                    'BoxStyle', 'outline', 'Symbol', 'k+')

            hold on

            
            scatter(group_indices+rand(size(group_indices))./4-0.125, sum(data_matrix(:,first_mode_idx:end),2), 60, cmap(group_indices,:), 'filled', 'o')
            xlabel('group');
            ylabel('complexity')
            box on; grid on
            % 
            % hold on
            % for gi = 1:length(gA.groups)
            %     data2plot = group_data(gi,:);
            % 
            %     plot(gi+rand(1,length(data2plot))./5-0.5/5, data2plot, ...
            %         '.', 'MarkerSize', 20, 'Color', cmap(gi,:));
            % 
            %     xlim([0 length(gA.groups)+1])
            % end

            set(gcf, 'Color', 'w')
            fontsize(gcf, 14, 'points')

            export_fig([gA.group_analysis_dir filesep 'LOCOEFA_complexity' num2str(t_idx, '%02d') '_' num2str(3) '_' num2str(5) '.png']);
            %%
            clf
            data_matrix = k(:,2:51);
            data_matrix = bsxfun(@rdivide, data_matrix, area_norm(:));
            norm = bsxfun(@rdivide, data_matrix, sum(data_matrix,2));
            norm_cum = cumsum(norm,2);
            
            last_mode = [];
            for i = 1:size(norm_cum,1)
                val = find(norm_cum(i,:) < 0.85, 1, 'last');
                if isempty(val); val = 1; end
                last_mode(i) = val;
            end
            
            group_data = nan(length(gA.groups), 100000);
            group_counter = ones(length(gA.groups), 1);
            for s_idx = 1:size(k, 1)
                group_id = results(s_idx).group_id;
                group_data(group_id, group_counter(group_id)) = last_mode(s_idx);
                group_counter(group_id) = group_counter(group_id)+1;
            end

            clf
            boxplot(group_data', gA.group_final_names,...
                    'ColorGroup', [], 'Colors', gA.param.cmap,...
                    'BoxStyle', 'outline', 'Symbol', 'k+')

            hold on

            scatter(group_indices+rand(size(group_indices))./2-0.25, last_mode, 60, cmap(group_indices,:), 'filled', 'o')
            title('required modes to explain 85% of the shape')
            xlabel('group')
            ylabel('modes')
            grid on; box on

            set(gcf, 'Color', 'w')
            fontsize(gcf, 14, 'points')

            export_fig([gA.group_analysis_dir filesep 'LOCOEFA_complexity_modes' num2str(t_idx, '%02d') '_' num2str(3) '_' num2str(5) '.png']);

            %%

        end

        %% Visualise individual samples
        % Plots anterior-posterior fluorescence profiles over time
        function plotAPprofile(gA, s_idx) %NOT WORKING AT THE MOMENT
            %% plot AP profile
            load([gA.group_analysis_dir filesep 'fluo_data.mat']); % loads fluo data variable

            S = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
            clf
            cmap = turbo(size(S,3));
            for g_idx = unique([fluo_data{s_idx}.g_idx])
                for t_idx = unique([fluo_data{s_idx}.t_idx])
                    for page = 1:length(fluo_data{s_idx})
                        if fluo_data{s_idx}(page).t_idx == t_idx && fluo_data{s_idx}(page).g_idx == g_idx 
                            subplot(2, max([fluo_data{s_idx}.g_idx]), g_idx)
                            title(['gastruloid id' num2str(g_idx)])
                            hold on
                            
                            plot([fluo_data{s_idx}(page).ch1_APprofile - fluo_data{s_idx}(page).ch1_Background],...
                                'color', cmap(t_idx,:))
                            
                            box on; grid on
                            if(g_idx)==1; ylabel('BRA');end
                            if page==length(fluo_data{s_idx})
                                
                                if(g_idx)==max([fluo_data{s_idx}.g_idx])
                                    colormap(cmap);
                                    cb = colorbar;
                                    cb.Ticks = [0, 0.5, 1];
                                    cb.TickLabels = {num2str(1), num2str(ceil(max([fluo_data{s_idx}.t_idx])/2)), num2str(max([fluo_data{s_idx}.t_idx]))};
                                    cb.Label.String = {'frame'};
                                end
                            end
        
                            subplot(2, max([fluo_data{s_idx}.g_idx]), g_idx+max([fluo_data{s_idx}.g_idx]))
                            hold on
                            plot([fluo_data{s_idx}(page).ch2_APprofile - fluo_data{s_idx}(page).ch2_Background],...
                                'color', cmap(t_idx,:))
                            
                            box on; grid on
                            xlabel('AP position')
                            if(g_idx)==1; ylabel('SOX2');end
                            if page==length(fluo_data{s_idx})
                                if(g_idx)==max([fluo_data{s_idx}.g_idx])
                                    colormap(cmap);
                                    cb = colorbar;
                                    cb.Ticks = [0, 0.5, 1];
                                    cb.TickLabels = {num2str(1), num2str(ceil(max([fluo_data{s_idx}.t_idx])/2)), num2str(max([fluo_data{s_idx}.t_idx]))};
                                    cb.Label.String = {'frame'};
                                end
                            end
                        end
                    end
                end
            end

            set(gca, 'Color', 'w')
            export_fig([gA.examples_dir filesep num2str(s_idx, '%03d') '_expression_profile.png'])
        end
        
        % Creates example visualizations with overlays and profile plots
        function showExample(gA, s_idx) %NOT WORKING AT THE MOMENT

                gA.examples_dir = [gA.main_dir filesep num2str(gA.current_idx, '%02d') '_examples'];
                mkdir(gA.examples_dir);

                images_path = gA.surf_dir;
                object_mask_path = gA.segm_dir;
                body_cell_mask_path = gA.semant_dir;
                data_dir = gA.properties_dir;
                save_folder = gA.examples_dir;
                
                %%
                stack_list = dir([images_path filesep '*.tif']);
                object_mask_list = dir([object_mask_path filesep '*.tif']);
                body_cell_mask_list = dir([body_cell_mask_path filesep '*.tif']);
                
                %%
                extend = @(x) x(:);
                
                % open stacks
                S = openTiffStack([stack_list(s_idx).folder filesep stack_list(s_idx).name]);
                O = openTiffStack([object_mask_list(s_idx).folder filesep object_mask_list(s_idx).name]);
                B = openTiffStack([body_cell_mask_list(s_idx).folder filesep body_cell_mask_list(s_idx).name]);
                
                % open json file with morgana results
                json_file = [data_dir filesep stack_list(s_idx).name '_morpho_params.json'];
                fid = fopen(json_file); 
                raw = fread(fid,inf); 
                str = char(raw'); 
                fclose(fid); 
                D = jsondecode(str);
                
                final_folder = [save_folder filesep stack_list(s_idx).name];
                mkdir(final_folder)
                for t_idx = 1:size(S,3)
                
                    clf
                    I = imcomplement(uint8(S(:,:,t_idx)));
                    props = regionprops(O(:,:,t_idx), 'Centroid');
                    
                    C = round(B(:,:,t_idx)/128);
                    
                    P = O(:,:, t_idx);
                    BW = P>0;
                    E = imdilate(BW, strel('disk', 10));
                    R = imerode(BW, strel('disk', 5));
                    
                    P = imdilate(P, strel('disk', 10));
                    P(double(E+R)~=1) = 0;
                    
                    I = labeloverlay(I, C, 'Transparency', 0.8, 'Colormap', parula(2));
                    I = labeloverlay(I, P, 'Transparency', 0.2, 'Colormap', lines(100));
                    
                    % determine number of labels within data set
                    g_inds = unique(P(P>0));
                    
                    theta = 0; % to rotate 90 counterclockwise
                    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
                    for g_idx = g_inds'
                        try
                        data_id = find(cellfun(@(x) contains(x, ['img_' num2str(t_idx, '%04d') '_cropped' num2str(g_idx,'%02d')]), {D.mask_file}));
                    
                        x = [(D(data_id).midline(:,[2,1])-D(data_id).centroid([2,1])')+props(g_idx).Centroid]';
                        j = mean(x,2);
                        y = R*(x-j)+j;
                    
                        line = extend(y)';
                        I = insertShape(I, 'line', line, 'LineWidth', 4, 'Color', 'black', 'Opacity', 0.1);
                        end
                    end
                
                    imwrite(I, [final_folder filesep 'img_' num2str(t_idx, '%04d') '.tif'])
                end

                plotAPprofile(gA, s_idx)
        end

        %% Diagnostic and Quality Control
        % Creates montage of all wells in 96-well plate format with channel overlays
        function plotSummary96WellPlate(gA, t_idx, scale)

            %% TO DO: super slow
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)

            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end

            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end

            height = size(S, 1);
            width = size(S, 2);
           
           second_elem = @(x) x{2};
           try %if names with the usual format "01_A01"
                well_names = cellfun(@(x) second_elem(split(x,'_')), gA.sample_names, 'UniformOutput',false);
           catch %if names with the format "A01"
               well_names = cellfun(@(x) x, gA.sample_names, 'UniformOutput',false);
           end
           
           row_ids = cellfun(@(x) strfind(string('A':'H'), x(1)), well_names);
           try
                col_ids = cellfun(@(x) str2double(x(2:3)), well_names);
           catch
                col_ids = cellfun(@(x) str2double(x(2)), well_names);
           end

           rows = 8;
           columns = 12;

           montage_size = [rows, columns];
           M = zeros(height*montage_size(1), width*montage_size(2), 3, 'uint8');
           M = imresize(M, scale);
        
            
            for s_idx = 1:length(gA.sample_names)
        
        %             S = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
                    C = openTiffFijiHyperStack([gA.channel_dir filesep gA.sample_names{s_idx} '.tif']);
        %             O = openTiffStack([gA.segm_dir filesep gA.sample_names{s_idx} '_mask.tif']);
        %             B = openTiffStack([gA.semant_dir filesep gA.sample_names{s_idx} '_mask.tif']);

                    h_i = row_ids(s_idx);
                    w_i = col_ids(s_idx);
        
                    channel_cmap = [1 1 1;... % gray
                                    1 0 0;... % red
                                    1 1 0;... % yellow
                                    0 1 1;... % cyan   
                                    0 1 0;... % green   
                                    0 0 1;... % blue
                                    1 0 1;... % magenta
                                    ];

                    channel_cmap = uint8(channel_cmap);

                    I = [];
                    for ch_idx = 1size(C,5):-1:1
                        ch = imadjust(uint8(C(:,:, 1, t_idx, ch_idx)));
                        % ch = uint8(C(:,:, 1, t_idx, ch_idx));
                        ch = cat(3, ch.*channel_cmap(ch_idx,1), ch.*channel_cmap(ch_idx,2), ch.*channel_cmap(ch_idx,3));
                        if isempty(I)
                            I = ch;
                        else
                            I = imfuse(ch, I, 'blend');
                        end
                    end

                    I = imresize(I, scale);
                    
                    M((h_i-1)*height*scale+1:(h_i-1)*height*scale+size(I,1), (w_i-1)*width*scale+1:(w_i-1)*width*scale+size(I,2),:) = I;
            end
        
            
            ML = M;
            for s_idx = 1:length(gA.sample_names)
                
                    h_i = row_ids(s_idx);
                    w_i = col_ids(s_idx);
                    
                    position = [(w_i-1)*width*scale+width*scale*0.01, (h_i-1)*height*scale+height*scale*0.01];
                   try %if names with the usual format "01_A01"
                        name = split(gA.sample_names{s_idx},'_');
                        name = [name{1} filesep name{2}];
                   catch %if names with the format "A01"
                        name = gA.sample_names{s_idx};
                   end
                    % w_i for column based experiments; h_i for row based
                    % experiments
                    ML = insertText(ML, position, name, 'FontSize', 50*scale,'AnchorPoint', 'LeftTop', 'BoxColor', round(cmap(w_i,:)*255));
            end
            imshow(ML)
            
        
        imwrite(ML, [save_dir filesep 'plate_summary_' num2str(t_idx, '%04d') '.png'])
        end

        % Shows segmentation masks overlaid on images in 96-well format
        function plotSummaryMasks96WellPlate(gA,  t_idx, scale)
            %%
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            groups = gA.groups;

            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end

            if isempty(gA.param.cmap)
                cmap = lines(length(groups));
            else
                cmap = gA.param.cmap;
            end

            height = size(S, 1);
            width = size(S, 2);
           
           second_elem = @(x) x{2};
           try %if names with the usual format "01_A01"
                well_names = cellfun(@(x) second_elem(split(x,'_')), gA.sample_names, 'UniformOutput',false);

           catch %if names with the format "A01"
               well_names = cellfun(@(x) x, gA.sample_names, 'UniformOutput',false);

           end
           
           row_ids = cellfun(@(x) strfind(string('A':'H'), x(1)), well_names);
           
           try
                col_ids = cellfun(@(x) str2double(x(2:3)), well_names);
           catch
                col_ids = cellfun(@(x) str2double(x(2)), well_names);
           end

           rows = 8;
           columns = 12;

           montage_size = [rows, columns];
           M = zeros(height*montage_size(1), width*montage_size(2), 3, 'uint8');
           M = imresize(M, scale);
            
            for s_idx = 1:length(gA.sample_names)
                    
                    S = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
%                     C = openTiffStack([gA.channel_dir filesep gA.sample_names{s_idx} '.tif']);
                    O = openTiffStack([gA.segm_dir filesep gA.sample_names{s_idx} '_mask.tif']);
        %             B = openTiffStack([gA.semant_dir filesep gA.sample_names{s_idx} '_mask.tif']);
                    
                    h_i = row_ids(s_idx);
                    w_i = col_ids(s_idx);
        
                    I = imadjust(uint8(S(:,:,t_idx)));
                    m = O(:,:,t_idx);
                    % w_i for column based experiments; h_i for row based experiments
                    L = labeloverlay(I, m>0, "Colormap", cmap(w_i,:));
                    if size(L, 3) == 1; L = cat(3, L, L, L); end

                    L = imresize(L, scale);
                    M((h_i-1)*height*scale+1:(h_i-1)*height*scale+size(L,1), (w_i-1)*width*scale+1:(w_i-1)*width*scale+size(L,2),:) = L;

            end
        
            
            ML = M;
            for s_idx = 1:length(gA.sample_names)
                
                    h_i = row_ids(s_idx);
                    w_i = col_ids(s_idx);
                    
                    position = [(w_i-1)*width*scale+width*scale*0.01, (h_i-1)*height*scale+height*scale*0.01];


                   try %if names with the usual format "01_A01"
                        name = split(gA.sample_names{s_idx},'_');
                        name = [name{1} filesep name{2}];
                   catch %if names with the format "A01"
                        name = gA.sample_names{s_idx};
                   end
                    % w_i for column based experiments; h_i for row based experiments
                    ML = insertText(ML, position, name, 'FontSize', 50*scale,'AnchorPoint', 'LeftTop', 'BoxColor', round(cmap(w_i,:)*255));
            end
            imshow(ML)
            
        
        imwrite(ML, [save_dir filesep 'mask_diagnosis_' num2str(t_idx, '%04d') '.png'])
        end

        % Creates montage of channel images for dish-based experiments
        function plotSummaryPictures(gA, t_idx, scale)
            
            %%
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)

            % use the size of the biggest image
            height = 0;
            width = 0;
            for idx = 1:length(gA.sample_names)
                clc
                disp(100*idx./length(gA.sample_names))
                S = openTiffStack([gA.surf_dir filesep gA.sample_names{idx} '.tif'], t_idx);

                height = max([height, size(S, 1)]);
                width = max([width, size(S, 2)]);
                
            end

            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end

            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end

%             height = size(S, 1);
%             width = size(S, 2);
           
           [rows, columns] = number_subplots(length(gA.sample_names));

           montage_size = [rows, columns];
           M = zeros(height*montage_size(1), width*montage_size(2), 3, 'uint8');
           M = imresize(M, scale);
            
            for s_idx = 1:length(gA.sample_names)
                    clc
                    disp(100*s_idx./length(gA.sample_names))
                    
                    [h_i, w_i] = ind2sub(montage_size, s_idx);
                    
                    if ~isempty(gA.channel_dir)
                        C = openTiffFijiHyperStack([gA.channel_dir filesep gA.sample_names{s_idx} '.tif'], 0, t_idx);
                    else
                        C = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif'], t_idx);
                    end
        
                    channel_cmap = [1 1 1;... % gray
                                    1 0 0;... % red
                                    1 1 0;... % yellow
                                    0 1 1;... % cyan   
                                    0 1 0;... % green   
                                    0 0 1;... % blue
                                    1 0 1;... % magenta
                                    ];

                    channel_cmap = uint8(channel_cmap);

                    I = [];
                    for ch_idx = size(C,4):-1:1
                        ch = imadjust(uint8(C(:, :, 1, ch_idx)));
                        ch = cat(3, ch.*channel_cmap(ch_idx,1), ch.*channel_cmap(ch_idx,2), ch.*channel_cmap(ch_idx,3));
                        if isempty(I)
                            I = ch;
                        else
                            I = imfuse(ch, I, 'blend');
                        end
                    end
                    
                    I = imresize(I, scale);

                    M(floor((h_i-1)*height*scale+1):floor((h_i-1)*height*scale+size(I,1)), floor((w_i-1)*width*scale+1):floor((w_i-1)*width*scale+size(I,2)),:) = I;

            end
        
            
            ML = M;
            for s_idx = 1:length(gA.sample_names)
                
                    [h_i, w_i] = ind2sub(montage_size, s_idx);
                    position = [(w_i-1)*width*scale+width*scale*0.01, (h_i-1)*height*scale+height*scale*0.01];
                    % name = split(gA.sample_names{s_idx},'_');
                    % name = [name{1} filesep name{2}];

                    % regex_res = cellfun(@(x) regexp(gA.sample_names{s_idx}, [x '\d+|' x '_'] )>0, gA.groups, 'UniformOutput', false); % determine group id
                    % group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last');
                    % name = gA.group_final_names{group_id};
                    regex_res = cellfun(@(x) contains(gA.sample_names{s_idx}, x), gA.groups);
                    group_id = find(regex_res, 1, 'last');
                    name = gA.group_final_names{group_id};

                    ML = insertText(ML, position, name, 'FontSize', 50*scale,'AnchorPoint', 'LeftTop', 'BoxColor', round(cmap(group_id,:)*255));
            end
            imshow(ML)
            
        
        imwrite(ML, [save_dir filesep 'dish_summary_' num2str(t_idx, '%04d') '.png'])

        end
        
        % Shows segmentation quality for dish-based experiments
        function plotSummaryMasks(gA, t_idx, scale)
            
            %%
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)
%             S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);

            % use the size of the biggest image
            height = 0;
            width = 0;
            for idx = 1:length(gA.sample_names)
                clc
                disp(100*idx./length(gA.sample_names))
                S = openTiffStack([gA.surf_dir filesep gA.sample_names{idx} '.tif'], t_idx);

                height = max([height, size(S, 1)]);
                width = max([width, size(S, 2)]);
                
            end


            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end

            if isempty(gA.param.cmap)
                cmap = lines(length(gA.groups));
            else
                cmap = gA.param.cmap;
            end

%             height = size(S, 1);
%             width = size(S, 2);
           
           [rows, columns] = number_subplots(length(gA.sample_names));

           montage_size = [rows, columns];
           M = zeros(height*montage_size(1), width*montage_size(2), 3, 'uint8');
           % M = 255*ones(height*montage_size(1), width*montage_size(2), 3, 'uint8');
           M = imresize(M, scale);
            
            for s_idx = 1:length(gA.sample_names)
                    clc
                    disp(100*s_idx./length(gA.sample_names))
                    
                    S = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif'], t_idx);
%                     C = openTiffStack([gA.channel_dir filesep gA.sample_names{s_idx} '.tif']);
                    O = openTiffStack([gA.segm_dir filesep 'raw_mask' filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '_cp_masks.tif']);
        %             B = openTiffStack([gA.semant_dir filesep gA.sample_names{s_idx} '_mask.tif']);
                    
                    [h_i, w_i] = ind2sub(montage_size, s_idx);

                    regex_res = cellfun(@(x) contains(gA.sample_names{s_idx}, x), gA.groups);
                    group_id = find(regex_res, 1, 'last');

                    % when we forget the word treatment
                    % regex_res = cellfun(@(x) regexp(gA.sample_names{s_idx}, ['\d+_' x '_\d+'] )>0, gA.groups, 'UniformOutput', false); % determine group id
                    % group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last');

                    I = imadjust(uint8(S));
                    % I = imcomplement(I);
                    m = O(:,:,1); % before tIdx
                    m = imgaussfilt(m, 5)>0.5;
                    m = imresize(m, size(I));

                    if ~isempty(gA.cellprob_segm_dir)

                        MC_t_idx = ones(size(I, 1), size(I, 2));
                        try 
                            % mask_filename = [gA.cellprob_segm_dir filesep 'raw_mask' filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '.tif_mask.tif'];
                            mask_filename = [gA.cellprob_segm_dir  filesep gA.sample_names{s_idx} '_mask.tif'];
                            MC_t_idx = imread(mask_filename);
                            MC_t_idx = imresize(MC_t_idx, [size(I, 1) size(I, 2)], 'nearest');
                        catch
                            disp([gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '.tif_mask.tif' ' not computed']);
                        end

                        c = imdilate(bwperim(MC_t_idx==1), strel('disk', 10));
                        m(c(:)>0) = 1;
                    end
                    
                    L = labeloverlay(I, m>0, "Colormap", cmap(group_id,:));
                    if size(L, 3) == 1; L = cat(3, L, L, L); end

                    L = imresize(L, scale);
                    M((h_i-1)*height*scale+1:(h_i-1)*height*scale+size(L,1), (w_i-1)*width*scale+1:(w_i-1)*width*scale+size(L,2),:) = L;

            end
        
            
            ML = M;
            for s_idx = 1:length(gA.sample_names)
                
                    [h_i, w_i] = ind2sub(montage_size, s_idx);
                    position = [(w_i-1)*width*scale+width*scale*0.01, (h_i-1)*height*scale+height*scale*0.01];
                    name = split(gA.sample_names{s_idx},'_');
                    name = name{2};
%                     name = [name{1} filesep name{3}];
% 
%                     regex_res = cellfun(@(x) regexp(gA.sample_names{s_idx}, [x '\d+|' x '_'] )>0, gA.groups, 'UniformOutput', false); % determine group id
%                     group_id = find(cellfun(@(x) ~isempty(x), regex_res), 1, 'last');
                    regex_res = cellfun(@(x) contains(gA.sample_names{s_idx}, x), gA.groups);
                    group_id = find(regex_res, 1, 'last');

                    ML = insertText(ML, position, name, 'FontSize', 50*scale,'AnchorPoint', 'LeftTop', 'BoxColor', round(cmap(group_id,:)*255));
            end
            imshow(ML)
            
        
        imwrite(ML, [save_dir filesep 'mask_diagnosis_' num2str(t_idx, '%04d') '.png'])

        end

        % Runs full diagnostic visualization suite for 96-well plates
        function plotDiagnosis96WellPlate(gA, scale) %make diagnosis functions for no plate mode
            if ~exist('scale', 'var') || isempty(scale)
                scale = 0.1;
            end
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);

            parfor t_idx = 1:size(S, 3)
                 gA.plotSummaryMasks96WellPlate(t_idx, scale)
            end

            parfor t_idx = 1:size(S, 3)
                 gA.plotSummary96WellPlate(t_idx, scale)
            end
        end
        
        % Runs full diagnostic visualization suite for dish experiments
        function plotDiagnosis(gA, scale) %make diagnosis functions for no plate mode
            if ~exist('scale', 'var') || isempty(scale)
                scale = 0.1;
            end
            S = openTiffStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
        
            parfor t_idx = 1:size(S, 3)
                 gA.plotSummaryMasks(t_idx, scale)
            end

            parfor t_idx = 1:size(S, 3)
                 gA.plotSummaryPictures(t_idx, scale)
            end
        end
    
        % Combines multiple analysis plots into single concatenated image
        function concatenateImagePlot(gA)
            concatimagelist={}; 

            v = openTiffFijiHyperStack([gA.surf_dir filesep gA.sample_names{1} '.tif']);
            for t_index = 1:size(v,4)
            
                numstring = num2str(t_index, '%02d');
                indir = [gA.group_analysis_dir filesep '*' numstring '*'];
                outdir = gA.main_dir;
                infiles = dir(indir);
                imagelist= {};
                
                for idx = 1:length(infiles)
                  I= imread([infiles(idx).folder filesep infiles(idx).name]);
                  imagelist{idx}=I;
                end
            
                if t_index==1
                    finalsize = [size(imagelist{1},1), size(imagelist{1},2)];
                end
                
                for idx= 1:length(infiles)
                    imagelist{idx}=imresize(imagelist{idx}, finalsize);
                end 
                
                combinedimages= cat(2, imagelist{:});
                
                concatimagelist{t_index} = combinedimages;
            end 
            
            combinedimagesintime = cat(1, concatimagelist{:});
            imwrite(combinedimagesintime, [outdir filesep 'combinedImage.png'])
    
        end 
    
        % Creates individual mask visualization for single sample
        function plotSingleMask(gA, s_idx, t_idx, colour, scale)
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)

            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end


            S = openTiffStack([gA.surf_dir filesep gA.sample_names{s_idx} '.tif']);
            O = openTiffStack([gA.segm_dir filesep 'raw_mask' filesep gA.sample_names{s_idx} '_t' num2str(t_idx, '%04d') '_cp_masks.tif']);          

            if ~exist('colour', 'var') || isempty(colour)
                regex_res = cellfun(@(x) contains(gA.sample_names{s_idx}, x), gA.groups);
                group_id = find(regex_res, 1, 'last');
                cmap = gA.param.cmap;
                colour = cmap(group_id,:);
            end

            I = imadjust(uint8(S(:,:,t_idx)));
            m = O(:,:,1); % before tIdx
            m = imgaussfilt(m, 5)>0.5;
            m = imresize(m, size(I));

            if ~isempty(gA.cellprob_segm_dir)
                c = openTiffStack([gA.cellprob_segm_dir filesep filesep gA.sample_names{s_idx} '_mask.tif']);
                c = imdilate(bwperim(c(:,:,1)), strel('disk', 10));
                m(c(:)>0) = 1;
            end
            
            L = labeloverlay(I, m>0, "Colormap", colour);
            if size(L, 3) == 1; L = cat(3, L, L, L); end

            L = imresize(L, scale);

            imwrite(L, [save_dir filesep gA.sample_names{s_idx} '_' num2str(t_idx, '%04d') '.png'])

        end

        % Creates individual channel visualization for single sample
        function plotSingleImage(gA, s_idx, t_idx, scale)
            save_dir = [gA.main_dir filesep 'mask_diagnosis'];
            mkdir(save_dir)

            if ~exist('t_idx', 'var') || isempty(t_idx); t_idx = 1; end

            if ~exist('scale', 'var') || isempty(scale)
                scale = 1;
            end


            C = openTiffFijiHyperStack([gA.channel_dir filesep gA.sample_names{s_idx} '.tif'], t_idx);
            channel_cmap = [1 1 1;... % gray
                            1 0 0;... % red
                            1 1 0;... % yellow
                            0 1 1;... % cyan   
                            0 1 0;... % green   
                            0 0 1;... % blue
                            1 0 1;... % magenta
                            ];

            channel_cmap = uint8(channel_cmap);

            I = [];
            for ch_idx = 1:-1:1%size(C,5):-1:1
                ch = imadjust(uint8(C(:, :, 1, ch_idx)));
                ch = cat(3, ch.*channel_cmap(ch_idx,1), ch.*channel_cmap(ch_idx,2), ch.*channel_cmap(ch_idx,3));
                if isempty(I)
                    I = ch;
                else
                    I = imfuse(ch, I, 'blend');
                end
            end
            
            I = imresize(I, scale);

            imwrite(I, [save_dir filesep gA.sample_names{s_idx} '_' num2str(t_idx, '%04d') '.png'])

        end
    end

end