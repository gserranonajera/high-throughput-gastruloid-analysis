function focus_Nikon_experiments_multistack(filename, outfolder, exp_name, workers, shift_idx, param)

    mkdir([outfolder])
    
    % Read stack sizes from info file.
    info = imfinfo(filename);
    info_text = info(1).ImageDescription;
    
    regex_template = {'(?<=','[^0-9]*)[0-9]*\.?[0-9]+'};

    firstElem = @(x)(str2double(x{1}));

    width = info(1).Width;
    height = info(1).Height;

    if contains(info_text, 'slices')
        depth = firstElem(regexp(info_text, [regex_template{1},'slices',regex_template{2}], 'match'));
    else
        depth = 1;
    end

    if contains(info_text, 'frames')
        num_frames = firstElem(regexp(info_text, [regex_template{1},'frames',regex_template{2}], 'match'));
    else
        num_frames = 1;
    end

    time_points = 1:num_frames;

    %% Gather parameters
    if ~exist('shift_idx', 'var') || isempty(shift_idx)
        shift_idx = 0;
    end

    if ~exist('param', 'var') || isempty(param)
        param = {
                         16,... % kernel_size
                         101,... % sq_side
                         5,... % overlaps
                         [0],... % z_averaging
                         1,... % channel
                         };
    end
    
    kernel_size = param{1};
    sq_side = param{2};
    overlaps = param{3};
    z_averaging = param{4};
    channel = param{5};
   
    
    %% start parallel pool
    pool = gcp('nocreate');
    if isempty(pool)
       parpool(workers);
    elseif pool.NumWorkers~=workers
       delete(pool);
       parpool(workers);
    end
    
    % split time points for the parallel workers for the images to be produced in order
    time_points_to_compute=cell(workers,1);
    for i_ind=1:workers
    %         time_points_to_compute{i_ind}=(i_ind+time_interval(1)-1):workers:time_interval(2);
          time_points_to_compute{i_ind}=time_points(i_ind:workers:num_frames);
    end
    
    img_vol_full = openTiffFijiHyperStack(filename);
    img_vol_full = squeeze(img_vol_full(:,:,:,:,channel));
    
    final_vol = zeros(height, width, num_frames, 'uint8');
    heightMapInTime = zeros(height, width, num_frames);
    img_list = {};
    hmap_list = {};
    %% Compute the heightmap

    if depth > 1
        
        % start work in all the workers
        parfor worker_ind=1:workers
        
            %% run surface finding algorithm independetly for each time point
            for t_ind=time_points_to_compute{worker_ind}
                
                % Read in raw input file.
%                 img_vol = tiffreadVolume(filename, 'PixelRegion', {[1 height], [1, width], [depth*(t_ind-1)+1 depth*t_ind]});
                img_vol = img_vol_full(:,:,:,t_ind);
                % Store the heigh map of each repetion
                heightMap_stack = zeros(height,width,overlaps);
                
                for disIdx = 1:overlaps % repetions
                
                    % every repetion dispace the columns. The mean of the maps
                    % is the final height map
                    displacement = (disIdx-1)*round(sq_side/overlaps)+1;
                    if disIdx>1 % avoid 0s in the final matrix
                        heightMap_stack(:,:,disIdx) = heightMap_stack(:,:,disIdx-1); 
                    end
                
                    % crop columns in the stack
                    for i = displacement:sq_side:size(img_vol,1)
                        for j = displacement:sq_side:size(img_vol,2)
                
                            final_i = min([i+sq_side-1, size(img_vol,1)]);
                            final_j = min([j+sq_side-1, size(img_vol,2)]);
                
                            % initilize vectors
                            flist = zeros(1,size(img_vol,3)); %focus
                            ilist = zeros(1,size(img_vol,3)); %intensity
                
                            % visit each slice in the column
                            for z = 1:size(img_vol,3)
                
                                subI = squeeze(img_vol(i:final_i, j:final_j, z));
                
                                % Squared gradient (Eskicioglu95)
                                Ix = diff(subI, 1, 2);
                                f = Ix.^2;
                                f = sum(f(:)); % measure the focus of the plane            
                                flist(z) = f; % measure the focus of the plane                    
                                ilist(z) = sum(subI(:)); % measure the average int. of the plane
                
                            end
                
                            % the location of the max. focus is the height
                            [~, h] = max(flist);
                            heightMap_stack(i:final_i, j:final_j,disIdx) = h;

                        end
                    end
                end
                
                % final heightmap
                heightMap = round(mean(heightMap_stack,3));
                kernel = ones(kernel_size)./kernel_size^2;
                heightMap=round(imfilter(heightMap,kernel,'same'));
                
                % Produce the surface. The image is the average of n planes in the surface
                vol = zeros(size(heightMap,1), size(heightMap,2), length(z_averaging), 'uint8');
                counter = 1;
                for z=z_averaging
                    for i = 1:height
                        for j = 1:width
                            zlevel = min([depth, max([1,heightMap(i,j)+z])]);
                            vol(i,j,counter) = img_vol(i,j,zlevel);
                        end
                    end
                    counter = counter + 1;
                end
    
                surface_00 = uint8(mean(vol,3));
%                 hmap_img = uint8(rescale(heightMap, 0, 255));
            
                img_list{worker_ind}{t_ind} = surface_00;
                hmap_list{worker_ind}{t_ind} = heightMap;

            end % t_ind

        end % workers

        for worker_ind=1:workers
            for t_ind = time_points_to_compute{worker_ind}
                final_vol(:,:,t_ind) = img_list{worker_ind}{t_ind};
                heightMapInTime(:,:,t_ind) = hmap_list{worker_ind}{t_ind};
            end
        end
        savetofile(heightMapInTime, [outfolder filesep exp_name '_hmap.mat'])
        saveTiffTimeStack(final_vol,[outfolder filesep exp_name '.tif'])

    else % if there is no z slices there is no need for surface finding

        heightMapInTime = ones(size(img_vol_full));
        savetofile(heightMapInTime, [outfolder filesep exp_name '_hmap.mat'])
        saveTiffTimeStack(img_vol_full,[outfolder filesep exp_name '.tif'])

%         % start work in all the workers
%         parfor worker_ind=1:workers
%         
%             for t_ind=time_points_to_compute{worker_ind}
%                 
%                 % Read in raw input file.
%                 surface_00 = tiffreadVolume(filename, 'PixelRegion', {[1 height], [1, width], [depth*(t_ind-1)+1 depth*t_ind]});
%                 heightMap = ones(size(surface_00));
% 
%                 img_list{worker_ind}{t_ind} = surface_00;
%                 hmap_list{worker_ind}{t_ind} = heightMap;
%             
%             end
% 
%         end
% 
%         for worker_ind=1:workers
%             for t_ind = time_points_to_compute{worker_ind}
%                 final_vol(:,:,t_ind) = img_list{worker_ind}{t_ind};
%                 heightMapInTime(:,:,t_ind) = hmap_list{worker_ind}{t_ind};
%             end
%         end
%         savetofile(heightMapInTime, [outfolder filesep exp_name '_hmap.mat'])
%         saveTiffTimeStack(final_vol,[outfolder filesep exp_name '.tif'])
    end

end % function