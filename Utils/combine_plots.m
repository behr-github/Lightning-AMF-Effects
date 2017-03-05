function combine_plots( varargin )
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
for a=1:numel(figs)
    nfigs = nfigs + sum(isgraphics(figs(a).Children, 'axes'));
    axinds = find(isgraphics(figs(a).Children,'axes'));
    inds = [repmat(a,numel(axinds),1), axinds(:)];
    figinds = cat(1, figinds, inds);
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
    
    f = figinds(a,1); c = figinds(a,2);
    copy_ax = copyobj(figs(f).Children(c), spfig);
    copy_ax.Position = pos;
end



end

