function [ rdiff ] = reldiff( A, B, varargin )
%RDIFF = RELDIFF( A, B ) Computes the relative difference of A - B.
%   Convenience function to make it require less typing to compute a
%   relative difference between A and B as (A - B)./B. Multiply the result
%   by 100 to convert to percent difference.
%
%   RDIFF = RELDIFF( A, B, TRUE ) will reshape the output into a column
%   vector.
%
%   RDIFF = RELDIFF( A, B, 'avg' ) will use the average of A and B in the
%   denominator rather than B. This is good for percent differences where
%   neither A nor B is "right" or "first" (i.e. you're not doing a percent
%   change or percent error).

E = JLLErrors;

forcevec = false;
denom_mode = 'first';

for a=1:numel(varargin)
    if (isnumeric(varargin{a}) || islogical(varargin{a})) && isscalar(varargin{a})
        forcevec = varargin{a};
    elseif ischar(varargin{a})
        denom_mode = varargin{a};
    else
        E.badinput('Optional arguments to RELDIFF must be either strings or scalar numbers/logicals.');
    end
end

if strcmpi(denom_mode, 'first')
    denom = B;
elseif strcmpi(denom_mode, 'avg')
    denom = 0.5 * (A+B);
else
    E.badinput('The only denominator modes allowed for RELDIFF are ''first'' and ''avg''');
end

rdiff = (A - B)./denom;

if forcevec
    rdiff = rdiff(:);
end

end

