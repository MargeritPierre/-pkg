function animation(frames,file,varargin) 
% EXPORT A SERIES OF FRAMES INTO A VIDEO
supportedOutputFormats = {'.avi' '.mp4' '.gif'} ;
if nargin<2 || isempty(file) 
    [file,path] = uiputfile(supportedOutputFormats,'SAVE THE FRAMES AS VIDEO OR GIF','animation.mp4') ;
    if path==0 ; return ; end
    file = [path filesep file] ;
end
file = string(file) ;
ext = regexp(file,supportedOutputFormats,'once') ;
if ~any(cat(1,ext{:})) ; file = string(file) + ".avi" ; end

% Convert to cell array of images
if isnumeric(frames) ; frames = num2cell(frames,[1 2 3]) ; end
if isstruct(frames) ; frames = {frames.cdata} ; end

% Build the animation file
[~,~,extension] = fileparts(file) ;
switch extension
    case ".gif" % GIF ANIMATION
    % Default options
        varargin = [...
                    {'LoopCount'}, {Inf} , ... infinite loop
                    {'DelayTime'}, {1/20} , ... frame rate
                    varargin];
    % Recording
        if size(frames{1},3)==1 % grayscale images
            imwrite(cat(4,frames{:}),file,"gif",varargin{:}) ;
        else
            for fr = 1:numel(frames)
                [A,map] = rgb2ind(frames{fr},256);
                if fr==1 
                    imwrite(A,map,file,"gif",varargin{:});
                else
                    imwrite(A,map,file,"gif","WriteMode","append",varargin{:});
                end
            end
        end
    otherwise % VIDEO ANIMATION
    % Default options
        switch extension
            case '.mp4'
                writerObj = VideoWriter(file,'MPEG-4') ;
                writerObj.Quality = 100 ;
                writerObj.FrameRate = 30 ;
            case '.avi'
                writerObj = VideoWriter(file,'Uncompressed AVI') ;
                writerObj.FrameRate = 30 ;
        end
    % User-specified options
        for arg = 1:2:numel(varargin)-1 
            writerObj.(varargin{arg}) = varargin{arg+1} ;
        end
    % Recording...
        open(writerObj) ;
        for fr = 1:numel(frames)
            writeVideo(writerObj,frames{fr}) ;
        end
        close(writerObj) ;
end

% Animation record


end