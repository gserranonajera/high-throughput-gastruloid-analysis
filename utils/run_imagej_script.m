function run_imagej_script(script_file, imagej_path, options, varargin)

    if ~exist('imagej_path','var') || isempty(imagej_path)
        if ispc
            %imagej_path = 'Y:\Users\Guillermo\software\fiji-win64_matlab\Fiji.app\ImageJ-win64.exe';
            imagej_path = 'C:\Users\gs714.BLUE\Desktop\Fiji.app\ImageJ-win64.exe';
        elseif isunix
            imagej_path='/home/gserranonajera/programs/Fiji.app/ImageJ-linux64';
        end
    end

    if ~exist('options','var') || isempty(options)
%         options='';
        options=' --no-splash -batch';
    end

    args=[];
    l=length(varargin);
    if l>0
        args='"';
        for i=1:length(varargin)
            args=[args varargin{i} ' '];
        end
        args=[args '"'];
    end

    if ispc
        args = strrep(args,'\','\\');
    end

    command=[imagej_path ' ' options ' ' script_file ' ' args];
    disp(command)
    system(command);
end