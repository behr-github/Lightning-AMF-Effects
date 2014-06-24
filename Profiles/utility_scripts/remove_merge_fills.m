function [ data, utc, alt, lon, lat, fills ] = remove_merge_fills( Merge, field )
%remove_merge_fills(Merge, field) Returns "field", utc, alt, lat, lon (properly
%scaled) with fill values in "field" replaced with NaNs.  Also returns the
%logical matrix of which entries were fills as the fifth output.

data = eval(sprintf('Merge.Data.%s.Values',field));
utc = Merge.Data.UTC.Values;
alt = Merge.Data.ALTP.Values;
lat = Merge.Data.LATITUDE.Values;
lon = Merge.Data.LONGITUDE.Values - 360;

fill_val = eval(sprintf('Merge.Data.%s.Fill',field));
ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;

fills = data == fill_val | data == ulod | data == llod;
data(fills) = NaN;
utc(fills) = NaN;
lat(fills) = NaN;
lon(fills) = NaN;
alt(fills) = NaN;

end

