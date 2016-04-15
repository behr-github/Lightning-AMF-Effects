function [ rdiff ] = reldiff( A, B, forcevec )
%RDIFF = RELDIFF( A, B ) Computes the relative difference of A - B.
%   Convenience function to make it require less typing to compute a
%   relative difference between A and B as (A - B)./B. Multiply the result
%   by 100 to convert to percent difference.
%
%   RDIFF = RELDIFF( A, B, TRUE ) will reshape the output into a column
%   vector.

if nargin < 3;
    forcevec = 0;
end

rdiff = (A - B)./B;

if forcevec
    rdiff = rdiff(:);
end

end

