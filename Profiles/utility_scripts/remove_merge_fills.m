%function [ data, utc, alt, lon, lat, fills ] = remove_merge_fills( Merge, field, varargin )
function [ data, varargout ] = remove_merge_fills( Merge, field, varargin )
%remove_merge_fills(Merge, field) Returns "field", utc, alt, lat, lon (properly
%scaled) with fill values in "field" replaced with NaNs.  Also returns the
%logical matrix of which entries were fills as the fifth output.  If you
%need to specify UTC, altitude, latitude, or longitude field names, do so
%as the parameters 'time', 'alt', 'lat', and 'lon'.

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('field',@isstr);
p.addParameter('time','UTC',@isstr);
p.addParameter('alt','ALTP',@isstr);
p.addParameter('lat','LATITUDE',@isstr);
p.addParameter('lon','LONGITUDE',@isstr);
p.addParameter('DEBUG_LEVEL',1,@(x) (isscalar(x) && isnumeric(x)));

p.parse(Merge,field,varargin{:});
pout = p.Results;
Merge = pout.Merge;
field = pout.field;
utcfield = pout.time;
altfield = pout.alt;
lonfield = pout.lon;
latfield = pout.lat;

E = JLLErrors;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

data = Merge.Data.(field).Values;
fill_val = Merge.Data.(field).Fill;
try
    ulod = Merge.metadata.upper_lod_flag;
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField')
        ulod = fill_val;
        if DEBUG_LEVEL > 0; fprintf('remove_merge_fills: ulod field does not exist\n'); end
    end
end

try
    llod = Merge.metadata.lower_lod_flag;
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField')
        llod = fill_val;
        if DEBUG_LEVEL > 0; fprintf('remove_merge_fills: llod field does not exist\n'); end
    end
end

if isnumeric(fill_val) && isscalar(fill_val)
    fills = data == fill_val | data == ulod | data == llod;
    data(fills) = NaN;
elseif strcmpi(fill_val,'n/a')
    % do nothing - this means that there should be absolutely no data
    % missing because no fill value was assigned in the Merge file. Usually
    % only happens for UTC.
else
    E.badvar('fill_val','Is not a numeric scalar value or "N/A"');
end
% Let the user specify longitude as the main input, in which case we need
% to correct the sign
if ~isempty(regexpi(field, 'longitude', 'ONCE'))
    data = lon_fix(data);
end
if nargout > 5; varargout{5} = fills; end

if nargout > 1; 
    utc = eval(sprintf('Merge.Data.%s.Values',utcfield)); 
    %utcfills = eval(sprintf('Merge.Data.%s.Fill',utcfield));
    %utc(utc==utcfills) = NaN; % There shouldn't ever be fill values in the UTC field
    varargout{1} = utc;
end
if nargout > 2; 
    alt = eval(sprintf('Merge.Data.%s.Values',altfield)); 
    altfills = eval(sprintf('Merge.Data.%s.Fill',altfield));
    alt(alt==altfills) = NaN;
    varargout{2} = alt;
end
if nargout > 4; 
    lat = eval(sprintf('Merge.Data.%s.Values',latfield)); 
    latfills = eval(sprintf('Merge.Data.%s.Fill',latfield));
    lat(lat==latfills) = NaN;
    varargout{4} = lat;
end
if nargout > 3; 
    lon = eval(sprintf('Merge.Data.%s.Values',lonfield)); 
    lonfills = eval(sprintf('Merge.Data.%s.Fill',lonfield));
    lon(lon==lonfills) = NaN;
    lon = lon_fix(lon);
    varargout{3} = lon;
end



end

function lon = lon_fix(lon)
lon(lon>180) = lon(lon>180) - 360;
end