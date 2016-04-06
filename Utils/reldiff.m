function [ rdiff ] = reldiff( A, B, forcevec )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3;
    forcevec = 0;
end

rdiff = (A - B)./B;

if forcevec
    rdiff = rdiff(:);
end

end

