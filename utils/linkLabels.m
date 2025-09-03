function linkLabels(path2label, out_directory, name, param, final_dim, link_labels)

    % links labels using basins in watershed!!

    if ~exist('param', 'var') || isempty(param) % these might need to change with different mask reduction sizes
        param = struct();
        param.filter_size = 11*4; %11
        param.close_size = 2*4; % 2
        param.gauss_std = 2*4; %2
        param.min_area = 250*4;
    end

    if ~exist('link_labels', 'var') || isempty(link_labels)
        link_labels = 1;
    end

    
%     label = openTiffStack(path2label);
    mask_files = dir(path2label);
    M = imread([mask_files(1).folder filesep mask_files(1).name]);
    label = zeros(size(M, 1), size(M, 2), length(mask_files));
    for t_idx = 1:length(mask_files)
        label(:,:,t_idx) = imread([mask_files(t_idx).folder filesep mask_files(t_idx).name]);
    end

    resized_label = zeros(final_dim);

    if link_labels % link labels based on basins
    
        new_label = double(label>0);
        final_label = zeros(size(new_label));
    
        basin_map = sum(new_label(:,:,:),3);
        ks = param.filter_size;
        kernel = ones(ks)./ks.^2;
        basin_map = imfilter(basin_map, kernel);
        basin_map = imgaussfilt(basin_map, param.gauss_std);
        basin_map = -basin_map;
        border_map = double(watershed(basin_map));
        
        for t_id = 1:size(label,3)
            new_label(:,:,t_id) = new_label(:,:,t_id).*border_map;
            ids = unique(extend(new_label(:,:,t_id)));
            ids(ids==0) = [];
            if ~isempty(ids)
                for i = ids'
                    BW = new_label(:,:,t_id)==i;
                    BW = bwareaopen(BW, param.min_area);
                    BW = double(imclose(BW, strel('disk', param.close_size)));
                    final_label(:,:,t_id) = final_label(:,:,t_id)+BW.*i;
                    resized_label(:,:,t_id) = imresize(final_label(:,:,t_id), final_dim(1:2), 'nearest');
                end
            end
        end

    else % don't link labels based on basins
        for t_id = 1:size(label,3)
            resized_label(:,:,t_id) = imresize(label(:,:,t_id), final_dim(1:2), 'nearest');
        end
    end

    %%
    vol = resized_label;
    file_name = [out_directory filesep name];

    fiji_descr = ['ImageJ=1.52p' newline ...
            'images=' num2str(1) newline... 
            'slices=' num2str(size(vol,3)) newline...
            'frames=' num2str(1) newline... 
            'hyperstack=false' newline...
            'mode=gray' newline...  
            'loop=false' newline...  
            'min=0.0' newline...      
            'max=256'];  % 8bit image
    
    t = Tiff(file_name,'w');
    tagstruct = {};
    tagstruct.ImageLength = size(vol,1);
    tagstruct.ImageWidth = size(vol,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % gray
    tagstruct.BitsPerSample = 8; % 8bit
    tagstruct.SamplesPerPixel = 1;% gray
    tagstruct.Compression = Tiff.Compression.LZW;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    tagstruct.ImageDescription = fiji_descr;
    for slice = 1:size(vol,3)
        t.setTag(tagstruct)
        t.write(uint8(vol(:,:,slice)));
        t.writeDirectory(); % saves a new page in the tiff file
    end
    
    t.close() 

end