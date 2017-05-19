function [  ] = label_axis_with_letter( label, varargin )
%LABEL_AXIS_WITH_LETTER Label a given axis with a string
%   LABEL_AXIS_WITH_LETTER( LABEL ) Places the string LABEL at the top left
%   corner of the current axes. If LABEL is a scalar number, it is
%   converted to a lowercase letter matching its place in the alphabet and
%   that letter, enclosed in parentheses, is used as the label.
%
%   Parameter arguments:
%       axis - give a handle to an axis to label instead of the current
%       one.
%
%       xshift - additional displacement from the axes in the x direction
%       as a fraction of the width of the plot. Default is 0.12.
%
%       yshift - additional displacement in the y direction as xshift is in
%       the x direction. Default is 0.
%
%       fontweight - any valid method of specifying a font weight to the
%       TEXT function. Default is 'b' (i.e. bold).
%
%       fontcolor - any valid color specification; changes the color of the
%       text. Default is 'k' (black).
%
%       fontsize - number specifying the font size in points. Default is
%       16.

E = JLLErrors;
if isnumeric(label) && isscalar(label)
    if label > 26
        warning('Numeric value of LABEL exceeds number of letters in alphabet');
    end
    label = sprintf('(%s)',char(label+96));
if ~ischar(label)
    E.badinput('LABEL must be a string or a scalar number')
end

p = inputParser;
p.addParameter('axis', gobjects(0));
p.addParameter('xshift', 0.12);
p.addParameter('yshift', 0);
p.addParameter('fontweight', 'b');
p.addParameter('fontcolor', 'k');
p.addParameter('fontsize', 16);

p.parse(varargin{:});
pout = p.Results;

ax = pout.axis;
shift_percent_x = pout.xshift;
shift_percent_y = pout.yshift;
fontweight = pout.fontweight;
fontcolor = pout.fontcolor;
fontsize = pout.fontsize;

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

text(x_left, y_top, label, 'color', fontcolor, 'fontweight', fontweight, 'fontsize', fontsize, 'parent', ax);

end

