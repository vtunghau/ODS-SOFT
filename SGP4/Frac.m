%--------------------------------------------------------------------------
% 
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function [res] = Frac(x)

res = x-floor(x);

