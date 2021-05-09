%% Satellite probagation
TLE2RV
clc
%close all
clear all
format long g
global AuxParam eopdata PC Cnm Snm const swdata

%% Load initial data

SAT_Const
load DE430Coeff.mat
load Epoch.mat
load State.mat
load TD.mat
load Kep.mat
load BStar.mat

PC = DE430Coeff;
Y0   = Y(1,:)';

%% Load File 
% read Earth gravity field coefficients
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
fclose(fid);

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data
fid = fopen('sw19571001.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);


%% Spacecraft probagation setup
Step   = 1;
Sample = (Mjd_Epoch(2)-Mjd_Epoch(1))*86400/Step;

%% AuxParam Setup
AuxParam.Mjd_UTC    = Mjd_Epoch(1);
AuxParam.n          = 30;
AuxParam.m          = 30;
AuxParam.sun        = 0;
AuxParam.moon       = 0;
AuxParam.planets    = 0;
AuxParam.sRad       = 0;
AuxParam.drag       = 0;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 0;
AuxParam.TD1        = TD.TD1(1);
AuxParam.TD2        = TD.TD2(1);
AuxParam.Bstar      = BStard(1);
AuxParam.Kep        = Kep(1,:);

Eph = Ephemeris(Y0,Sample,Step);
Y2  = State(const.GM_Earth,Kep(2,:),-((Mjd_Epoch(2)-Mjd_Epoch(1))*86400 - Eph(length(Eph),1)));

plot3(Eph(:,2),Eph(:,3),Eph(:,4));              hold on
plot3(Y2(1),Y2(2),Y2(3),'or');                  hold on


%% Update State
Y_S         = State(const.GM_Earth,Kep(1,:),Eph(length(Eph),1));
Kep3        = Element(const.GM_Earth,Eph(length(Eph),2:7))
Kep2        = Element(const.GM_Earth,Y_S);
KepChange   = Kep3 - Kep2;
KepF        = Kep3;

E           = EccAnom(Kep2(6), Kep2(2));
v           = 2*atan(sqrt((1+Kep2(2))/(1-Kep2(2)))*tan(E/2));
v           = 2*pi + v - KepChange(5) ;
KepF(6)     = mod(v2M(v,Kep2(2)),2*pi) - sin(Kep2(3))*KepChange(4);
Y_F         = State(const.GM_Earth,KepF,0);



Y2  = State(const.GM_Earth,Kep(2,:),-((Mjd_Epoch(2)-Mjd_Epoch(1))*86400 - Eph(length(Eph),1)));

Error           = ((Y_F' - Y2').^2)';
Error_dis       = sqrt(sum(Error(1:3,:)))'
Error_dis_dt    = Error_dis/((Mjd_Epoch(2)-Mjd_Epoch(1))*86400);
Error_vel       = sqrt(sum(Error(4:6,:)))';

plot3(Y_F(1),Y_F(2),Y_F(3),'xk');               hold on
plot3(Y2(1),Y2(2),Y2(3),'og');                  hold on


%% Update more frequently



