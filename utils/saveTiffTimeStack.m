function saveTiffTimeStack(stack, file_name)

    fiji_descr = ['ImageJ=1.52p' newline ...
                'images=' num2str(size(stack, 3)) newline... 
                'slices=' num2str(1) newline...
                'frames=' num2str(size(stack,3)) newline... 
                'hyperstack=false' newline...
                'mode=8bit' newline...  
                'loop=false' newline...  
                'min=0.0' newline...      
                'max=256'];  % 8bit image
    
    t = Tiff(file_name,'w');
    tagstruct = {};
    tagstruct.ImageLength = size(stack,1);
    tagstruct.ImageWidth = size(stack,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % rgb
    tagstruct.BitsPerSample = 8; % 8bit
    tagstruct.SamplesPerPixel = 1;% rgb
    tagstruct.Compression = Tiff.Compression.LZW;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    tagstruct.ImageDescription = fiji_descr;
    for frame = 1:size(stack,3)
        t.setTag(tagstruct)
        t.write(im2uint8(stack(:,:,frame)));
        t.writeDirectory(); % saves a new page in the tiff file
    end
    t.close() 

end