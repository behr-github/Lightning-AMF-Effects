function varargout = combine_plots( varargin )
%COMBINE_PLOTS Combine individual plots into a single figure
%   COMBINE_PLOTS() Combined all open figures into a single figure
%   keeping the number of axes along each dimension as close to
%   equal as possible.
%
%   COMBINE_PLOTS(FIGS) Combine only the figures specified by the
%   array of handles FIGS.
%
%   COMBINE_PLOTS( ___, 'dims', [y x] ) Combined with either previous
%   syntax, specify the dimensions as vertical x horizontal (as in 
%   subplot)

E = JLLErrors;
p = inputParser;

p.addOptional('figs',[]);
p.addParameter('dims',[]);
p.parse(varargin{:});
pout = p.Results;
figs = pout.figs;

if isempty(figs)
    figs = get(0,'children');
elseif any(~isgraphics(figs,'figure'))
    E.badinput('FIGS must contain only figure handles');
end

% Calculate how many axes we'll need in the subplot
nfigs = 0;
figinds = [];
n_obj_types = 3; % how many object types we'll try to copy. Right now 3: axes, legend, and colorbar
inds = nan(1,n_obj_types+1);
for a=1:numel(figs)
    nfigs = nfigs + sum(isgraphics(figs(a).Children, 'axes'));
    inds(1) = a;
    % Go through the list of children for this figure. Assume that any
    % associated objects (e.g. colorbar) are listed before the associated
    % axis, so find each type of object we want to copy, then when we find
    % an axis instance, reset.
    for b=1:numel(figs(a).Children)
        if isgraphics(figs(a).Children(b), 'colorbar')
            inds(3) = b;
        elseif isgraphics(figs(a).Children(b), 'legend')
            inds(4) = b;
        elseif isgraphics(figs(a).Children(b), 'axes')
            inds(2) = b;
            figinds = cat(1, figinds, inds);
            inds = nan(1,n_obj_types+1);
            inds(1) = a;
        end
    end
end

% Make the subplots roughly even in each dimension, preferring width over
% height
if isempty(pout.dims)
    x = ceil(sqrt(nfigs));
    y = ceil(nfigs/x);
else
    if ~isnumeric(pout.dims) || numel(pout.dims) > 2 || any(pout.dims < 1)  || any(mod(pout.dims,1) > 0)
        E.badinput('''dims'' must be a two-element positive vector of whole numbers')
    elseif prod(pout.dims) < nfigs
        E.badinput('''dims'' contains insufficient plots (%d axes to combine found)',nfigs)
    end
    y = pout.dims(1);
    x = pout.dims(2);
end

spfig = figure;
for a=1:nfigs
    sp = subplot(y,x,a);
    sp.Units = 'Normalized';
    pos = sp.Position;
    delete(sp);
    
    f = figinds(a,1); 
    cc = ~isnan(figinds(a,:)); cc(1) = false; % never copy the whole figure (which is the handle index in the first column)
    c = figinds(a,cc);
    copy_ax = copyobj(figs(f).Children(c), spfig);
    new_ax = isgraphics(copy_ax, 'axes');
    copy_ax(new_ax).Position = pos;
end

if nargout > 0
    varargout{1} = spfig;
end


end

