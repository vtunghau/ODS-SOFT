%--------------------------------------------------------------------------
%
% Sun: Computes the Sun's geocentric position using a low precision
%      analytical series
%
% Input:
%   Mjd_TT    Terrestrial Time (Modified Julian Date)
% 
% Output:     
%   rSun      Solar position vector [m] with respect to the 
%             mean equator and equinox of J2000 (EME2000, ICRF)
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function rSun = Sun(Mjd_TT)

SAT_Const
pi2 = 2*pi;
ep  = const.Rad*(84381.412/3600);     % Obliquity of J2000 ecliptic
T   = (Mjd_TT-const.MJD_J2000)/36525; % Julian cent. since J2000

% Mean anomaly, ecliptic longitude and radius
M = pi2 * Frac( 0.9931267 + 99.9973583*T);                 % [rad]
L = pi2 * Frac( 0.7859444 + M/pi2 + ...
              (6892.0*sin(M)+72.0*sin(2.0*M)) / 1296.0e3); % [rad]
r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2*M);         % [m]

% Equatorial position vector
rSun = R_x(-ep) * [r*cos(L), r*sin(L), 0]';

