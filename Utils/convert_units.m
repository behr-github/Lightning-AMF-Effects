function [ data ] = convert_units( data, unit_in, unit_out )
%CONVERT_UNITS Convert data between defined units
%   Takes in a matrix or vector of data and two strings defining the unit
%   the data is in and the unit you want the data to be in.  Will convert
%   between the units as long as they are defined and in the same category.
%   You can also pass 'ls', 'list', or 'listunits' as the only input to see
%   what units it knows and their categories.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 June 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add units to a category (or make a new category here). The code will
% require that both the in and out unit be in the same category.
UnitCategories.mixing_ratios.units = {'ppt','pptv','ppb','ppbv','ppm','ppmv','ppp','pppv'};
UnitCategories.mixing_ratios.fxn = @mixing_ratio_conv;
UnitCategories.pressure.units = {'Pa', 'hPa', 'kPa'};
UnitCategories.pressure.fxn = @pressure_conv;

if ischar(data) && ismember(data, {'ls','list','listunits'})
    list_units(UnitCategories)
    data = 1;
    return
elseif ~isnumeric(data) 
    E.badinput('data is expected to be numeric')
elseif ~ischar(unit_in) || ~ischar(unit_out)
    E.badinput('The units are expected to be given as strings')
end



units_given = {unit_in, unit_out};
fns = fieldnames(UnitCategories);
both_same_cat = false;
for a = 1:numel(fns)
    in_cat = ismember(units_given, UnitCategories.(fns{a}).units);
    if all(in_cat)
        both_same_cat = true;
    elseif any(in_cat) && ~all(in_cat)
        in_cat = find_cat(UnitCategories, unit_in);
        out_cat = find_cat(UnitCategories, unit_out);
        E.badinput('The given units are not in the same category (%s in %s, %s in %s',unit_in,in_cat,unit_out,out_cat)
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONVERTING %%%%%
%%%%%%%%%%%%%%%%%%%%%%

unit_cat = find_cat(UnitCategories, unit_in);
conv_fxn = UnitCategories.(unit_cat).fxn;

in_conv = conv_fxn(unit_in);
out_conv = conv_fxn(unit_out);

data = data / in_conv * out_conv;

end

function list_units(UnitCategories)
fns = fieldnames(UnitCategories);
for a=1:numel(fns)
    units = strjoin(UnitCategories.(fns{a}).units, ', ');
    fprintf('%s: %s\n', fns{a}, units);
end
end

function cat_name = find_cat(UnitCategories, unit)
fns = fieldnames(UnitCategories);
for a=1:numel(fns)
    if ismember(unit, UnitCategories.(fns{a}).units);
        cat_name = fns{a};
        return
    end
end
end

function conv = mixing_ratio_conv(unit)
E=JLLErrors;
switch unit
    case 'ppp'
        conv = 1;
    case 'pppv'
        conv = 1;
    case 'ppm' 
        conv = 1e6;
    case 'ppmv'
        conv = 1e6;
    case 'ppb'
        conv = 1e9;
    case 'ppbv'
        conv = 1e9;
    case 'ppt'
        conv = 1e12;
    case 'pptv'
        conv = 1e12;
    otherwise
        E.badinput('Unit %s not recognized as a mixing ratio')
end
end

function conv = pressure_conv(unit)
E=JLLErrors;
switch unit
    case 'Pa'
        conv = 1;
    case 'hPa'
        conv = 1e-2;
    case 'kPa' 
        conv = 1e-3;
    otherwise
        E.badinput('Unit %s not recognized as a mixing ratio')
end
end