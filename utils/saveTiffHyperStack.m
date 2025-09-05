function saveTiffHyperStack(vol, file_name)
    
    %dimensions: [x, y, z, t, c] (3 because it is RGB)
    
    fiji_descr = ['ImageJ=1.52p' newline ...
                'images=' num2str(size(vol,3)*size(vol,4)*size(vol,5)) newline... 
                'slices=' num2str(size(vol,3)) newline...
                'frames=' num2str(size(vol,4)) newline...
                'channels=' num2str(size(vol,5)) newline...
                'hyperstack=true' newline...
                'mode=rgb' newline...  
                'loop=false' newline...  
                'min=0.0' newline...      
                'max=256'];  % 8bit image
    
    t = Tiff(file_name,'w8'); %w8 to write bigtiffs
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
    for slice = 1:size(vol,3)
        for frame = 1:size(vol,4)
            for channel = 1:size(vol,5)
                t.setTag(tagstruct)
                t.write(uint8(vol(:,:,slice,frame,channel)));
                t.writeDirectory(); % saves a new page in the tiff file
            end
        end
    end
    t.close() 
end