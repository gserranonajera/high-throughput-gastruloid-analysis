function tiffStack=openTiffStack(fileName, page_idx)
    TifFile=fileName;
    InfoImage=imfinfo(TifFile);
    width=InfoImage(1).Width;
    height=InfoImage(1).Height;
    depth=length(InfoImage);
    
    
    if ~exist('page_idx', 'var') || isempty(page_idx)
        tiffStack=zeros(height,width,depth,'double');
        for i=1:depth
           tiffStack(:,:,i)=imread(TifFile,'Index',i,'Info',InfoImage);
        end
    else
        tiffStack=zeros(height,width,1,'double');
        tiffStack = imread(TifFile,'Index',page_idx,'Info',InfoImage);
    end
end