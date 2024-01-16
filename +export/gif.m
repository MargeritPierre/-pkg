function gif(frames,file,varargin) 
% EXPORT A SERIES OF FRAMES INTO A GIF
if nargin<2 || isempty(file) 
    [file,path] = uiputfile('frames.gif','SAVE THE FRAMES INTO A GIF') ;
    if path==0 ; return ; end
    file = [path filesep file] ;
end
file = string(file) ;
if ~any(regexp(file,".gif",'once')) ; file = string(file) + ".gif" ; end

if iscell(frames) ; frames = cat(4,frames{:}) ; end

% Default options
varargin = [...
            {'LoopCount'}, {Inf} , ... infinite loop
            {'DelayTime'}, {1/20} , ... frame rate
            varargin];

if size(frames,3)==1 % grayscale images
    imwrite(frames,file,"gif",varargin{:}) ;
else
    for fr = 1:size(frames,4)
        [A,map] = rgb2ind(frames(:,:,:,fr),256);
        if fr==1 
            imwrite(A,map,file,"gif",varargin{:});
        else
            imwrite(A,map,file,"gif","WriteMode","append",varargin{:});
        end
    end
end


end