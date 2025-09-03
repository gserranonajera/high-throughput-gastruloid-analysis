function run_morgana_analysis_flourescence(image_stack_dir, mask_stack_dir, save_folder, plate_mode)
%{
Function to call morgana. 
1) Each stack (images and masks) needs to be decomposed into individual images. 
2) Call morgana pipeline in python (run_morgana_pipeline.py)
3) Move the json files with the results to the main results folder
4) Try to delete all the unnecessary data

Important: if the label has disconnected pieces (or there are very small
labels) morgana crashes
%}
    image_stack_list = dir([image_stack_dir filesep '*.tif']);
    mask_stack_list = dir([mask_stack_dir filesep '*.tif']);
    morgana_script = 'Y:\Users\Guillermo\software\repository\gastruloid_analysis_functions\run_morgana_pipeline.py';
  
    

    current_folder = [save_folder];
    current_folder_images = [current_folder filesep 'images'];
    current_folder_masks = [current_folder filesep 'masks'];

    mkdir(current_folder_images)
    mkdir(current_folder_masks)
    parfor s_idx = 1:length(image_stack_list)
        clc
        disp([num2str(s_idx) ' ' image_stack_list(s_idx).name])

        final_file = [save_folder filesep, ...
            image_stack_list(s_idx).name '_morpho_params.json'];

        if exist(final_file, 'file')
            continue
        end

        S = openTiffFijiHyperStack([image_stack_list(s_idx).folder filesep image_stack_list(s_idx).name], []);
        M = openTiffStack([mask_stack_list(s_idx).folder filesep mask_stack_list(s_idx).name]);
    
        for tIdx = 1:size(S,4)
            I = uint8(squeeze(S(:,:,:,tIdx,:)));
%             m = uint8(M(:,:,tIdx));
            
            % smooth signal to avoid issues with morgana
            m = imgaussfilt(double(M(:,:,tIdx)), 15)>0.5;
            m = bwareaopen(m, 1000, 4);
            m = uint8(m);
            m = imclose(m, strel('disk', 25)); % remove corners
            % plate modes assumes there is one gastruloid, the biggest
            % object
            if plate_mode 
                props = regionprops(m, 'area');
                if ~isempty(props)
                    [~, max_id] = max([props.Area]);
                    m = uint8(m==max_id)*max_id;
                end            
            else
                bw = zeros(size(m), 'uint8');
                for label = unique(m(m>0))'
                    % remove small objects and smooth the labels for morgana
%                     bw = bw+imopen(bwareaopen(m==label, 30), strel('disk', 30));
                    bw = bw + m==label;%imgaussfilt(double(m==label), 20)>0.5;
                end
                
                m = m.*uint8(bw);
            end

            %imwrite(I, [current_folder_images filesep image_stack_list(s_idx).name '_' num2str(tIdx, '%04d') '.tif']);
            saveTiffStack(I, [current_folder_images filesep image_stack_list(s_idx).name '_' num2str(tIdx, '%04d') '.tif'])
            imwrite(m, [current_folder_masks  filesep image_stack_list(s_idx).name '_' num2str(tIdx, '%04d') '.tif']);

        end
    end

    command = ['conda activate morgana &',...
       'python ' morgana_script, ...
       ' ' current_folder_images, ...
       ' ' current_folder_masks, ...
       ' ' 'aggregated_results', ...
       ];
    
    system(command)

    res_file_a = [current_folder_images filesep,...
                'splitObjects' filesep 'result_segmentation' filesep,...
                'splitObjects_morpho_params.json'];

    res_file_b = [current_folder_images filesep,...
                'splitObjects' filesep 'result_segmentation' filesep,...
                'aggregated_results' '_fluo_intensity.json'];

%         res_file_b = [current_folder_images filesep,...
%             'splitObjects' filesep 'result_segmentation' filesep,...
%             image_stack_list(idx).name '_morpho_straight_params.json'];

    try
        copyfile(res_file_a, [save_folder filesep 'aggregated_results' '_morpho_params.json']);
        copyfile(res_file_b, [save_folder filesep 'aggregated_results' '_fluo_intensity.json']);
        rmdir(current_folder_images, 's');
        rmdir(current_folder_masks, 's');
    catch E
        disp(E)
    end
end
