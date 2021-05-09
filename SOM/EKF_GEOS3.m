%--------------------------------------------------------------------------
%
%  Initial Orbit Determination using Gauss and Extended Kalman Filter methods
%
% References:
%   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
%   Applications", Springer Verlag, Heidelberg, 2000.
%   
%   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
%   4th Edition, 20 13.
%
%   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
%
% Last modified:   2020/03/16   Meysam Mahooti
%--------------------------------------------------------------------------
clc
clear
format long g

addpath('C:\Users\VTBB3\Desktop\ODS\SOM\');
global const Cnm Snm AuxParam eopdata n_eqn PC

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end

% Model parameters
AuxParam = struct ('Mjd_UTC',0,'n',0,'m',0);

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);



load ..\Data\AZ.mat

AZ  = AZ(AZ(:,2)>-17,:);
Eph = AZ(:,6:11);
AZ  = AZ(:,[[1,4,5] 12]);
%% read observations

nobs = length(AZ);
obs = zeros(nobs,4);
for i=1:nobs
%     tline = fgetl(fid);
%     if ~ischar(tline)
%         break;
%     end
%     Y = str2num(tline(1:4));
%     M = str2num(tline(6:7));
%     D = str2num(tline(9:10));
%     h = str2num(tline(13:14));
%     m = str2num(tline(16:17));
%     s = str2num(tline(19:24));
%     az = str2num(tline(26:33));
%     el = str2num(tline(36:42));
%     Dist = str2num(tline(45:54));
    az= AZ(i,2);
    el= AZ(i,3);
    obs(i,1) = AZ(i,1);
    obs(i,2) = az;
    obs(i,3) = el;
    obs(i,4) = AZ(i,4);
    %obs(i,4) = 1e3*Dist;
end
% 
% fclose(fid);
sigma_az    = 0.01*const.Rad; % [rad]
sigma_el    = 0.01*const.Rad; % [rad]
sigma_range = 1;
load ..\Data\Inertial.mat

Rs          = Inertial.R_Sta';
LT          = Inertial.E;

Mjd1 = obs(1,1);
Mjd2 = obs(9,1);
Mjd3 = obs(18,1);


% [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
%                   Mjd1,Mjd2,Mjd3,Rs,Rs,Rs)
[r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
                     Mjd1,Mjd2,Mjd3,Rs,Rs,Rs)


Y0_apr = [r2;v2];
Y0_apr - Eph(9,:)';

Mjd0 = obs(1,1);

Mjd_UTC = obs(1,1);

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 71;
AuxParam.m       = 71;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides=0;
AuxParam.sRad = 0;
AuxParam.drag = 0;
AuxParam.Relativity = 0;
n_eqn  = 6;

Y = DEInteg(@Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
Eph(1,:)' - Y
P = zeros(6);
  
for i=1:3
    P(i,i)=1e8;
end
for i=4:6
    P(i,i)=1e3;
end



yPhi = zeros(42,1);
Phi  = zeros(6);

% Measurement loop
t = 0;

for i=1:nobs    
    % Previous step
    i
    t_old = t;
    Y_old = Y;
    
    % Time increment and propagation
    Mjd_UTC = obs(i,1);                       % Modified Julian Date
    t       = (Mjd_UTC-Mjd0)*86400.0;         % Time since epoch [s]
    
    [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
    [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;
        
    for ii=1:6
        yPhi(ii) = Y_old(ii);
        for j=1:6  
            if (ii==j) 
                yPhi(6*j+ii) = 1; 
            else
                yPhi(6*j+ii) = 0;
            end
        end
    end
    
    yPhi = DEInteg (@VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
    
    % Extract state transition matrices
    for j=1:6
        Phi(:,j) = yPhi(6*j+1:6*j+6);
    end
    Y = DEInteg (@Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
    
    % Topocentric coordinates
    theta = gmst(Mjd_UT1);                    % Earth rotation
    U = R_z(theta);
    r = Y(1:3);
    s = LT*(U*r-Rs);                          % Topocentric position [m]
    
    % Time update
    P = TimeUpdate(P, Phi);
        
    % Azimuth and partials
    [Azim, ~, dAds, ~] = AzElPa(s);     % Azimuth, Elevation
    dAdY = [dAds*LT*U,zeros(1,3)];
    
    % Measurement update
    [~, Y, P] = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
    
    % Elevation and partials
    r = Y(1:3);
    s = LT*(U*r-Rs);                          % Topocentric position [m]
    [Azim, Elev, dAds, dEds] = AzElPa(s);     % Azimuth, Elevation
    dEdY = [dEds*LT*U,zeros(1,3)];
    
    % Measurement update
    [~, Y, P] = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
    
%     %Range and partials
%     r = Y(1:3);
%     s = LT*(U*r-Rs);                          % Topocentric position [m]
%     Dist = norm(s); dDds = (s/Dist)';         % Range
%     dDdY = [dDds*LT*U,zeros(1,3)];
%     
%     % Measurement update
%     [K, Y, P] = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(i,1),'l');
    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT  = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Y0      = DEInteg (@Accel,0,-(obs(i,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);
    Eph(1,:)' - Y0
    Diff(i,:) = Eph(1,:)' - Y0;
end

rmpath('C:\Users\VTBB3\Desktop\ODS\SOM\');