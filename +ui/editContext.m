function H = editContext(obj,varargin)
%EDITCONTEXT create context menu for in-situ property modification

% Input check
    if nargin==0 ; obj = gcf ; end
    if ~isa(obj,'matlab.graphics.Graphics') ; error('The object should be a graphic handle.') ; end
    if mod(numel(varargin),2)~=0 ; error('Incorrect Nam-Value argument pair.') ; end

% Initialize the handles
    H = gobjects(0) ;

% Is there a need to create a context menu ?
    if isa(obj,'matlab.ui.container.ContextMenu') ; return ; end
    if isa(obj,'matlab.ui.container.Menu') ; return ; end
    if ~isprop(obj,'UIContextMenu') ; return ; end
    %disp(class(obj))
    
% Apply to children
    if isprop(obj,'Children') && ~isempty(obj.Children)
        % Apply (recusively) the function to children
        hh = arrayfun(@(child)pkg.ui.editContext(child,varargin{:}),obj.Children,'UniformOutput',false) ;
        H = [H hh{:}] ;
    end
    
% Get the settable properties & infos
    propInfos = set(obj) ;
    if isempty(propInfos) ; return ; end
    propNames = fieldnames(propInfos) ;
    if isempty(propNames) ; return ; end
    
% Check if the property can be simply set by menus
    values = cellfun(@obj.get,propNames,'UniformOutput',false) ;
    settable = cellfun(@ischar,values) ;
    settable = settable | cellfun(@isnumeric,values) ;
    propNames = propNames(settable) ;
    if isempty(propNames) ; return ; end
    
% Initialize the context menu
    fig = ancestor(obj,'figure') ;
    ui = uicontextmenu(fig) ;
    
% Create the menu tree
    sm = createMenuTree(propNames,ui) ;
    
% Add menus associated to the properties
    addmenus(obj,propNames,sm) ;

% Add the menu to the handles
    obj.UIContextMenu = ui ;
    H(end+1) = ui ;

end

function str = splitWords(str)
% Split a 'WordWithCapitals' to a {'Word' 'With' 'Capitals'}
    isCapital = [1 find(ismember(str,'A':'Z')) numel(str)+1] ;
    lengths = diff(isCapital) ;
    str = mat2cell(str,1,lengths) ;
    str(lengths==0) = [] ;
end

function str = value2str(value)
% Convert a value to a string
    nMax = 5 ;
    precision = 4 ;
    if numel(value)<=nMax
        str = mat2str(value,precision) ;
    else
        str = regexprep(mat2str(value(1:nMax),precision),']',' ... ]') ;
    end
end

function setProp(obj,propName)
% Open a dialog to set an object's property
    prompt = {'Modify the property value:'} ;
    title = [obj.Type '.' propName ] ;
    definput = {mat2str(obj.(propName))} ;
    dims = [10 50] ;
    resizeable = 'on' ;
    answer = inputdlg(prompt,title,dims,definput,resizeable) ;
    if isempty(answer) ; return ; end
    newValue = str2num(answer{1}) ;
    if isempty(newValue) ; newValue = answer{1} ; end
    obj.(propName) = newValue ;
end

function sm = createMenuTree(propNames,ui)
% Create the menu tree corresponding to the list of property names
    sm = gobjects(numel(propNames),1) ;
    % Build the tree of property names
        propWords = cell(numel(propNames),100) ;
        [propWords{:}] = deal('') ;
        for pp = 1:numel(propNames)
            words = splitWords(propNames{pp}) ;
            propWords(pp,1:numel(words)) = words ;
        end
    % Reduce the tree size
        emptyWords = cellfun(@isempty,propWords) ;
        propWords = propWords(:,~all(emptyWords,1)) ;
    % Reduce to unique branches
        for ww = 1:size(propWords,2)-1
            fullWords = join(propWords(:,1:ww),'',2) ;
            [~,ia,iu] = unique(fullWords) ;
            isUnique = sum(ia(iu)==1:size(propWords,1),1)==1 ;
            propWords(isUnique,ww) = join(propWords(isUnique,ww:end),'',2) ;
            [propWords{isUnique,ww+1:end}] = deal('') ; 
        end
    % Add all parent menus
        sm(:) = ui ; % initially target to the top tree branch
        for ww = 1:size(propWords,2)
            for pp = 1:size(propWords,1)
                if emptyWords(pp,ww) ; continue ; end
                fullWord = strjoin(propWords(pp,1:ww),'') ;
                [~,ind] = ismember(fullWord,get(sm,'Tag')) ;
                if ind==0 
                    sm(pp) = uimenu(sm(pp),'Text',propWords{pp,ww},'Tag',fullWord) ;
                else
                    sm(pp) = sm(ind) ;
                end
            end
        end
end

function addmenus(obj,propNames,sm)
    for pp = 1:numel(propNames)
        % Values
            possibleValues = set(obj,propNames{pp}) ; % Possible values
            currentValue = get(obj,propNames{pp}) ; % Current Value
        % Create the submenu
            if isempty(possibleValues) % Numeric data or any string
                ssm = uimenu(sm(pp),'Text',value2str(currentValue)) ;
                ssm.Position = 1 ;
            else % options
                ssm = gobjects(0) ;
                for vv = 1:numel(possibleValues)
                    val = possibleValues{vv} ;
                    switch class(val)
                        case 'char'
                            ssm(end+1) = uimenu(sm(pp),'Text',val) ;
                            if strcmp(currentValue,val) ; ssm(end).Checked = 'on' ; end
                        otherwise
                    end
                end
            end
        % Set the submenu callbacks
            sm(pp).UserData = obj ;
            sm(pp).Callback = @parentMenuCallback ;
            set(ssm,'UserData',obj) ; 
            set(ssm,'Tag',propNames{pp}) ; 
            set(ssm,'Callback',@menuCallback) ;
    end
end


function parentMenuCallback(src,evt)
    obj = src.UserData ;
    propName = src.Tag ;
    value = obj.(propName) ;
    if isempty(set(obj,propName)) % numeric value or any string
        ssm = src.Children(ismember({src.Children.Tag},propName)) ;
        ssm.Text = value2str(value) ;
    else % options
        [src.Children.Checked] = deal('off') ;
        set(src.Children(ismember({src.Children.Text},value)),'Checked','on') ;
    end
end


function menuCallback(src,evt)
    obj = src.UserData ;
    propName = src.Tag ;
    if isempty(set(obj,propName)) % numeric value or any string
        setProp(obj,propName) ;
    else % options
        [src.Parent.Children.Checked] = deal('off') ;
        src.Checked = 'on' ;
        obj.(propName) = src.Text ;
    end
end







