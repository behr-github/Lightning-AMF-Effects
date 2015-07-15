function [ varargout ] = campaign_wide_ops( campaign_name, species, op, varargin )
%CAMPAIGN_WIDE_OPS Operations performed on a campaign's worth of data.
%   When it is desirable to look at data for an entire campaign, it is very
%   inconvinient that each day is its own Merge file. This function can run
%   several different operations on all the data in the campaign. Takes at
%   least 3 arguments:
%       campaign_name - the name of the campaign to operate on. Must be
%       recognized by merge_field_names.
%
%       species - Essentially what field to operate on in the Merge file.
%       This will try first to find this in the Names structure returned by
%       merge_field_names (so no2_lif will get you the appropriate NO2 LIF
%       field for the given campaign), then try to find that field in the
%       Merge files for the campaign.
%
%       op - the operation to be carried out. These will be described in
%       more detail below. Some operations may require additional arguments
%       to be passed; see their descriptions below for this information.
%
%   This function will return one or more outputs, again this depends on
%   the operation.
%
%   Possible operations are:
%
%       'cat' - simply concatenate the designated species into one long
%       vector. Will return this vector, plus a vector of UTC times and a
%       cell array of dates describing the date/time of each measurement.
%       Requires no additional arguments.
%
%       'bin' - bins the data by GPS altitude. Requires 1 additional
%       argument, the bin width in km. Outputs the bin median values, bin
%       midpoint altitudes, and upper and lower quartiles for each bin.
%
%       'bin_rolling' - bins the data with rolling averaging windows (using
%       bin_rolling_vertical_profile). Need the bin width and bin spacing
%       (in km) passed as additional arguments. Outputs the bin median
%       values, bin midpoint altitudes, and upper and lower quartiles for
%       each bin.
%
%       'bin_pres' - bins by pressure, into OMI standard pressure bins.
%       Requires no additional arguments. Returns the bin median values,
%       bin pressure levels, and upper and lower quartiles.

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(3,Inf);

% Confirm that the operation requested is allowed, then check that enough
% arguments were passed.

allowed_ops =   {'cat','bin','bin_rolling','bin_pres'};
req_args =      [0,    1,    2,            0];

op = lower(op);
if ~ismember(op, allowed_ops)
    E.badinput('op %s is not one of the expected values: %s', op, strjoin(allowed_ops, ', '));
end
xx = strcmp(op, allowed_ops);
if numel(varargin) < req_args(xx)
    E.badinput('op %s requires %d additional arguments', op, req_args(xx));
end

% We'll check that species is a valid input later, we need a Merge file
% loaded to do that

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

[Names, ~, merge_dir] = merge_field_names(campaign_name);

F = dir(fullfile(merge_dir,'*.mat'));
Ms(numel(F)) = struct('Merge',[]);
for a=1:numel(F)
    if DEBUG_LEVEL > 0; fprintf('Loading %s\n', F(a).name); end
    load(fullfile(merge_dir, F(a).name),'Merge'); % brings the variable Merge into the workspace
    % Now we'll handle checking species against Names and Merge fields. On
    % successive loops, ensure that the field name is still defined for the
    % new Merge.
    if a == 1
        if isfield(Names, species)
            merge_field = Names.(species);
        elseif isfield(Merge.Data, species)
            merge_field = species;
        else
            E.badinput('species %s is not a defined field in Names or Merge.Data for the campaign %s', species, campaign_name);
        end
    elseif ~isfield(Merge.Data, merge_field)
        E.callError('inconsistent_merge','Later Merge.Data does not have the field %s', merge_field);
    end
    
    Ms(a).Merge = Merge;
end

% At this point, the operation specified must be called.
switch op
    case 'cat'
        [varargout{1}, varargout{2}, varargout{3}] = concatenate_merges(Ms, merge_field);
    case 'bin'
        [varargout{1}, varargout{2}, varargout{3}] = bin_merges(Ms, merge_field, varargin{1}, Names);
    case 'bin_rolling'
        [varargout{1}, varargout{2}, varargout{3}] = bin_rolling_merges(Ms, merge_field, varargin{1:2}, Names);
    case 'bin_pres'
        [varargout{1}, varargout{2}, varargout{3}] = bin_pressure_merges(Ms, merge_field, Names);
    otherwise
        E.badinput('Operation %s not recognized',op);
end


end

function [catted_species, utcs, dates] = concatenate_merges(Ms, merge_field)
catted_species = [];
utcs = [];
dates = {};
% Loop through all the merges, appending the requested species, times, and
% dates
for a=1:numel(Ms)
    catted_species = cat(2, catted_species, remove_merge_fills(Ms(a).Merge, merge_field));
    utcs = cat(2, utcs, Ms(a).Merge.Data.UTC.Values);
    n = numel(Ms(a).Merge.Data.UTC.Values);
    this_dates = repmat({Ms(a).Merge.metadata.date},1,n);
    dates = cat(2, dates, this_dates);
end
end

function [bin_vals, bin_alts, bin_quarts] = bin_merges(Ms, merge_field, bin_width, Names)
all_species = concatenate_merges(Ms, merge_field);
all_alts = concatenate_merges(Ms, Names.gps_alt);
[bin_vals, bin_alts, bin_quarts] = bin_vertical_profile(all_alts, all_species, bin_width);
end

function [bin_vals, bin_alts, bin_quarts] = bin_rolling_merges(Ms, merge_field, bin_width, bin_spacing, Names)
all_species = concatenate_merges(Ms, merge_field);
all_alts = concatenate_merges(Ms, Names.gps_alt);
[bin_vals, bin_alts, bin_quarts] = bin_rolling_vertical_profile(all_alts, all_species, bin_width, bin_spacing);
end

function [bin_vals, bin_pres, bin_quarts] = bin_pressure_merges(Ms, merge_field, Names)
all_species = concatenate_merges(Ms, merge_field);
all_pres = concatenate_merges(Ms, Names.pressure);
[bin_vals, bin_pres, bin_quarts] = bin_omisp_pressure(all_pres, all_species);
end