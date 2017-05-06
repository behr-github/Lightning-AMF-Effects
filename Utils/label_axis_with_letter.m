function [  ] = label_axis_with_letter( label, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;
if ~ischar(label)
    E.badinput('LABEL must be a string')
end

p = inputParser;
p.addParameter('axis', gobjects(0));
p.addParameter('xshift', 0.12);
p.addParameter('yshift', 0);
p.addParameter('fontweight', 'b');
p.addParameter('fontcolor', 'k');
p.addParameter('fontsize', 16);
p.addParameter('parent', gca);

p.parse(varargin{:});
pout = p.Results;

ax = pout.axis;
shift_percent_x = pout.xshift;
shift_percent_y = pout.yshift;
fontweight = pout.fontweight;
fontcolor = pout.fontcolor;
fontsize = pout.fontsize;
parent = pout.parent;

if ~isempty(ax)
    if ~ishandle(ax) || ~strcmp(get(ax,'type'), 'axes')
        E.badinput('AX (if given) must be a handle to axes')
    end
else
    ax = gca;
end

% Determine the top-left axis limits
xl = get(ax, 'xlim');
if strcmp(get(ax, 'xdir'), 'normal')
    x_left = min(xl);
else
    x_left = max(xl);
end
x_left = x_left - shift_percent_x * diff(xl);

yl = get(ax, 'ylim');
if strcmp(get(ax, 'ydir'), 'normal')
    y_top = max(yl);
else
    y_top = min(yl);
end
y_top = y_top + shift_percent_y * diff(yl);

text(x_left, y_top, label, 'color', fontcolor, 'fontweight', fontweight, 'fontsize', fontsize, 'parent', parent);

end

