function [  ] = label_subfigs( fig, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xx = find(strcmpi(varargin, 'ax'));
if ~isempty(xx)
    varargin(xx:xx+1) = [];
end

xx = isgraphics(fig.Children, 'axes');
figax = fig.Children(xx);
pos = cat(1, figax.Position);
[~,ord_ind] = sortrows(pos,[-2 1]); % sort by vertical position in reverse order, then by horizontal position in forward order

for a=1:numel(ord_ind)
    % add the subfigure letters. 
    label_axis_with_letter(a, 'ax', figax(ord_ind(a)), varargin{:});
end

end

