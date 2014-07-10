function [ varargout ] = state_outlines( varargin )
%state_outlines(states): Draws state outlines using MATLAB shapefiles
%   Plots outlines of states using MATLAB shape files, which makes them
%   easier to plot more stuff on top of than using m_map. Pass specific
%   state abbreviations, or 'all' (or no arguments) to draw to full US

not_states = false;
if nargin == 0;
    states = 'all';
elseif strcmpi(varargin{1},'not')
    states = varargin(2:end);
    not_states = true;
else
    states = varargin;
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
fnum = figure;
if nargout > 0; varargout{1} = fnum; end
if strcmpi(states,'all');
    for a=1:numel(usa)
        plot(usa(a).X, usa(a).Y,'color','k');
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
            plot(usa(b).X, usa(b).Y, 'color','k')
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
        plot(usa(xx).X, usa(xx).Y,'color','k');
        hold on
    end
end
end

