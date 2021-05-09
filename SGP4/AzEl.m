%--------------------------------------------------------------------------
%
% AzEl: Computes azimuth and elevation from local tangent coordinates
%
% Input:
%   s   Topocentric local tangent coordinates (East-North-Zenith frame)
% 
% Outputs:
%   A   Azimuth   [rad]
%   E   Elevation [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [A, E] = AzEl(s)

A = atan2(s(1),s(2));

if (A<0) 
    A = A+2*pi;
end

E = atan( s(3)/sqrt(s(1)^2+s(2)^2) );

