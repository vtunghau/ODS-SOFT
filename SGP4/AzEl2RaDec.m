function [RA DEC] = AzEl2RaDec(Az,El,lat,lon,Mjd_UT1)
% Programed by Darin C. Koblick 6/16/2009
%--------------------------------------------------------------------------
% Updated:                                                      Date:
% - quadrant check bug fix                                  1/22/2010
% - vectorized for speed                                    1/22/2010
%--------------------------------------------------------------------------
% External Function Call Sequence:
% [RA DEC] = AzEl2RaDec(0,0,0,-104,'1992/08/20 12:14:00')
%
% Worked Example pg. 262 Vallado
% [RA DEC] = AzEl2RaDec(210.8250667,23.8595052,39.007,-104.883,'1994/05/14 13:11:20.59856')
% [alpha_t,delta_t] = AzEl2RaDec(Beta,el,phi,lamda,'yyyy/mm/dd hh:mm:ss')
%
% Function Description:
%--------------------------------------------------------------------------
% AzEl2RaDec will take the Azimuth and Elevation in the local horizon
% reference frame, site latitude and longitude as well as a time in GMT
% and output the Right Ascension and Declination in the topocentric coordinate frame.
%
% Inputs:                                                       Format:
%--------------------------------------------------------------------------
% Local Azimuth Angle   (degrees)                               [N x 1]
% Local Elevation Angle (degrees)                               [N x 1]
% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))            [N x 1]
% Lon (Site Longitude in degrees -180:180 W(-) E(+))            [N x 1]
% UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)          [N x 1]
%
% Outputs:                                                      Format:
%--------------------------------------------------------------------------
% Topocentric Right Ascension (Degrees)   [N x 1]
% Topocentric Declination Angle (Degrees)                       [N x 1]
%
%
% External Source References:
% Fundamentals of Astrodynamics and Applications 
% D. Vallado, Second Edition
% Example 3-5. Finding Local Siderial Time (pg. 192) 
% Algorithm 28: AzElToRaDec (pg. 259)
% -------------------------------------------------------------------------
%Example 3-5
ThetaGMST = gmst(Mjd_UT1)*180/pi;
ThetaLST = ThetaGMST + lon;
%Algorithm 28
DEC = asind(sind(El).*sind(lat)+cosd(El).*cosd(lat).*cosd(Az));
LHA = atan2(-sind(Az).*cosd(El)./cosd(DEC), ...
    (sind(El)-sind(DEC).*sind(lat))./(cosd(DEC).*cosd(lat))).*(180/pi);
RA = mod(ThetaLST-LHA,360);