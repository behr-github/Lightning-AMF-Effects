function [ varargout ] = state_outlines( varargin )
%state_outlines(states): Draws state outlines using MATLAB shapefiles
%   Plots outlines of states using MATLAB shape files, which makes them
%   easier to plot more stuff on top of than using m_map. Pass specific
%   state abbreviations, or 'all' (or no arguments) to draw to full US
%
%   Pass a figure number as the first argument to draw the outlines on the
%   specified figure.  gcf returns a number, and so works just fine as this
%   argument.
%
%   Most color specifications can be passed as the second argument (or the
%   first if not also passing a figure number).  This will recolor the
%   state outlines from their default black.  Note that only 1 x 3 vectors
%   or 1 character strings will work.
%
%   Josh Laughner <joshlaugh5@gmail.com> Jul 2014

argin = varargin;

if isnumeric(argin{1}) && numel(argin{1}) == 1;
    fignum = argin{1};
    argin(1) = [];
elseif ishandle(argin{1}) && strcmp(get(argin{1},'Type'),'figure')
    fignum = get(argin{1},'Number');
    argin(1) = [];
else 
    fignum = 0;
end

if ischar(argin{1}) && length(argin{1}) == 1;
    colspec = argin{1}; 
    argin(1) = [];
elseif isnumeric(argin{1}) && length(argin(1))==3;
    colspec = argin{1};
    argin(1) = [];
else
    colspec = 'k';
end
    

not_states = false;
if nargin == 0 || isempty(argin);
    states = 'all';
elseif strcmpi(argin{1},'not')
    states = argin(2:end);
    not_states = true;
else
    states = argin;
end

state_abbrev = {'al','ak','az','ar','ca','co','ct','de','fl','ga','hi','id','il','in','ia','ks','ky','la','me','md','ma','mi','mn','ms','mo','mt','ne','nv','nh','nj','nm','ny','nc','nd','oh','ok','or','pa','ri','sc','sd','tn','tx','ut','vt','va','wa','wv','wi','wy'};
state_names = {'alabama','alaska','arizona','arkansas','california','colorado','connecticut',...
    'delaware','florida','georgia','hawaii','idaho','illinois','indiana','iowa','kansas',...
    'kentucky','louisiana','maine','maryland','massachusetts','michigan','minnesota','mississippi',...
    'missouri','montana','nebraska','nevada','new_hampshire','new_jersey','new_mexico','new_york','north_carolina','north_dakota',...
    'ohio','oklahoma','oregon','pennsylvania','rhode_island','south_carolina','south_dakota','tennessee','texas','utah','vermont','virginia','washington','west_virginia',...
    'wisconsin','wyoming'};

% Read the states shapefile and plot
usa = shaperead('usastatehi.shp');

if fignum < 1; fnum = figure;
else figure(fignum); fnum = fignum;
end

% Figure out is the axes are a map or not. This may break in matlab
% versions before 2014b, when the graphics object was introduced.
ch = get(figure(fignum),'children');
im_a_map = false;
for a=1:numel(ch)
    if strcmpi(get(ch(a),'Type'),'axes')
        im_a_map = ismap(ch(a));
        break
    end
end

if nargout > 0; varargout{1} = fnum; end
if strcmpi(states,'all');
    for a=1:numel(usa)
        if im_a_map
            linem(usa(a).Y, usa(a).X, 'color', colspec);
        else
            line(usa(a).X, usa(a).Y,'color',colspec);
        end
        hold on
    end
elseif not_states
    states_not_to_draw = cell(1,numel(states));
    for c=1:numel(states) 
        if length(states{c}) == 2
            xx = find(strcmpi(states{c},state_abbrev));
        else
            xx = find(strcmpi(states{c},state_names));
        end
        states_not_to_draw{c} = state_names{xx};
    end
    for b = 1:numel(state_names)
        if any(strcmpi(usa(b).Name,states_not_to_draw))
            continue
        else
            if im_a_map
                linem(usa(b).Y, usa(b).X, 'color', colspec);
            else
                line(usa(b).X, usa(b).Y, 'color',colspec)
            end
            hold on
        end
    end
    hold off
else
    for b = 1:numel(states);
        if length(states{b}) == 2
            xx = find(strcmpi(states{b},state_abbrev));
        else
            xx = find(strcmpi(states{b},state_names));
        end
        if im_a_map
            linem(usa(xx).Y, usa(xx).X, 'color', colspec);
        else
            line(usa(xx).X, usa(xx).Y,'color',colspec);
        end
        hold on
    end
end
end

