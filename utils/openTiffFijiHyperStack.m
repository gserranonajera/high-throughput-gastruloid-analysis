function image = openTiffFijiHyperStack(img_path, big_tiff, t_idx)
    
    % t_idx only works with big_tiff=False

    if  ~exist('big_tiff','var') || isempty(big_tiff)
        big_tiff = false;
    end

    if  ~exist('t_idx','var') || isempty(t_idx)
        t_idx = [];
    end

    info = imfinfo(img_path);
    info_text = info(1).ImageDescription;
    
    regex_template = {'(?<=','[^0-9]*)[0-9]*\.?[0-9]+'};

    firstElem = @(x)(str2double(x{1}));

    width = info(1).Width;
    height = info(1).Height;

    try 
        depth = firstElem(regexp(info_text, [regex_template{1},'slices',regex_template{2}], 'match'));
    catch
        depth = 1;
    end

    try
        num_images = firstElem(regexp(info_text, [regex_template{1},'images',regex_template{2}], 'match'));
    catch
        num_images = 1;
    end

    try
        num_channels = firstElem(regexp(info_text, [regex_template{1},'channels',regex_template{2}], 'match'));
    catch
        num_channels = 1;
    end

    try
        num_frames = firstElem(regexp(info_text, [regex_template{1},'frames',regex_template{2}], 'match'));
    catch
        num_frames = 1;
    end

    % sometimes necessary to create an info file per image
    if length(info) ~= num_images
        for i = 1:num_images
            info(i) = info(1);
        end
    end
    
    if big_tiff
        clc
        disp(['Loading big tiff...']);
        vol = tiffreadVolume(img_path);

        % faster than preallocating
        image = reshape(vol, [height, width, depth, num_frames, num_channels]);

        % somehow the command before does not work properly, need to
        % reorganise the stack with this loop
        counter = 1;
        % order ch -> z -> t
        for t = 1:num_frames
            for z = 1:depth
                for ch = 1:num_channels
                    clc
                    disp(['Loading TIFF stack: ' num2str(counter*100./num_images) '%'])
                    image(:, :, z, t, ch) = vol(:,:,counter);                   
                    counter = counter + 1;
                end
            end
        end

    else %no big_tiff
        
        probe = imread(img_path, 1, 'Info', info);
        
        if isempty(t_idx)
            disp(['Prealocating memory...'])
            if isa(probe, 'uint16')
                image = uint16(zeros(height, width, depth, num_frames, num_channels));
            
            elseif isa(probe, 'uint8')
                image = uint8(zeros(height, width, depth, num_frames, num_channels));
            end
    
            counter = 1;
            % order ch -> z -> t
            for t = 1:num_frames
                for z = 1:depth
                    for ch = 1:num_channels
                        clc
                        disp(['Loading TIFF stack: ' num2str(counter*100./num_images) '%'])
                        image(:, :, z, t, ch) = imread(img_path, counter, 'Info', info);                   
                        counter = counter + 1;
                    end
                end
            end

        else % with t_idx

            if isa(probe, 'uint16')
                image = uint16(zeros(height, width, depth, num_channels)); %only 1 frame
            
            elseif isa(probe, 'uint8')
                image = uint8(zeros(height, width, depth, num_channels)); %only 1 frame
            end
    
            %jump to the current timepoint
            counter = num_channels*depth*(t_idx-1)+1;
            % order ch -> z -> t
            for z = 1:depth
                for ch = 1:num_channels
                    image(:, :, z, ch) = imread(img_path, 'Index', counter, 'Info', info);                   
                    counter = counter + 1;
                end
            end
            


        end %t_idx

    end % big tiff
end