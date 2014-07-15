function [ data, utc, alt, lon, lat, fills ] = remove_merge_fills( Merge, field, varargin )
%remove_merge_fills(Merge, field) Returns "field", utc, alt, lat, lon (properly
%scaled) with fill values in "field" replaced with NaNs.  Also returns the
%logical matrix of which entries were fills as the fifth output.  If you
%need to specify UTC, altitude, latitude, or longitude field names, do so
%as the parameters 'time', 'alt', 'lat', and 'lon'.

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('field',@isstr);
p.addParamValue('time','UTC',@isstr);
p.addParamValue('alt','ALTP',@isstr);
p.addParamValue('lat','LATITUDE',@isstr);
p.addParamValue('lon','LONGITUDE',@isstr);

p.parse(Merge,field,varargin{:});
pout = p.Results;
Merge = pout.Merge;
field = pout.field;
utcfield = pout.time;
altfield = pout.alt;
lonfield = pout.lon;
latfield = pout.lat;

data = eval(sprintf('Merge.Data.%s.Values',field));
utc = eval(sprintf('Merge.Data.%s.Values',utcfield));
alt = eval(sprintf('Merge.Data.%s.Values',altfield));
lat = eval(sprintf('Merge.Data.%s.Values',latfield));
lon = eval(sprintf('Merge.Data.%s.Values',lonfield));
lon(lon>180) = lon(lon>180) - 360;

fill_val = eval(sprintf('Merge.Data.%s.Fill',field));
ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;

fills = data == fill_val | data == ulod | data == llod;
data(fills) = NaN;

end

