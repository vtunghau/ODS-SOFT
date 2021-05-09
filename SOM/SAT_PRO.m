%--------------------------------------------------------------------------
%
%   Exercise 2-4: Topocentric satellite motion
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
TLE2RV
clc
close all
clear all
format long g
global counter AuxParam eopdata PC Cnm Snm const swdata

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

load Epoch.mat
load State.mat
load TD.mat
load Kep.mat
load BStar.mat
Y0   = Y(1,:)';

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

% Spacecraft orbit


Step   = 1;
Sample = (Mjd_Epoch(2)-Mjd_Epoch(1))*86400/Step;


% Orbit
for Sec= 1:(Sample+1)
    
    Mjd_UTC = Mjd_Epoch(1) + (Sec-1)/(86400/Step); % Time
    dt = (Mjd_UTC-Mjd_Epoch(1))*86400;    % Time since epoch [s] 
    
    Y_S(Sec,:) = State(const.GM_Earth, Kep(1,:), (Sec-1)*Step);          % Inertial vector
end
counter = 0;
AuxParam.Mjd_UTC    = Mjd_Epoch(1);
AuxParam.n          = 180;
AuxParam.m          = 180;
AuxParam.sun        = 1;
AuxParam.moon       = 1;
AuxParam.planets    = 1;
AuxParam.sRad       = 0;
AuxParam.drag       = 1;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;
AuxParam.TD1        = TD.TD1(1);
AuxParam.TD2        = TD.TD2(1);
AuxParam.Bstar      = BStard(1);
AuxParam.Kep        = Kep(1,:);



Eph = Ephemeris(Y0,Sample,Step);
Error = ((Eph(:,2:7) - Y_S).^2)';
%plot(sqrt(sum(Error(4:6,:)))')
figure
%plot(Eph(:,2))
%hold on
%plot(Y_S(:,1))
%legend('Pos','Vel');
%Y2 = State(const.GM_Earth, Kep(2,:),-0.473313733935*60)
AuxParam.Mjd_UTC = Mjd_Epoch(2);
AuxParam.Kep = Kep(2,:);
Y2 = Ephemeris(Y(2,:)',60,-1);
% plot(Y2(:,1)/60 + (Mjd_Epoch(2)-Mjd_Epoch(1))*86400/60+1,Y2(:,2))

plot3(Eph(:,2),Eph(:,3),Eph(:,4)); hold on
plot3(Y2(:,2),Y2(:,3),Y2(:,4));
% Error = ((Eph(length(Eph),2:7) - Y2(1,2:7)).^2)';
% sqrt(sum(Error(1:3,:)))'
% ans/((Mjd_Epoch(2)-Mjd_Epoch(1))*86400)
% sqrt(sum(Error(4:6,:)))'



Kep3 = Element(const.GM_Earth,Eph(length(Eph),2:7));
%Kep3 = Element(const.GM_Earth,Y2(1,2:7));
Kep2 = Element(const.GM_Earth,Y_S(length(Y_S),:));
KepChange = Kep3 - Kep2;
KepF = Kep3;
E           = EccAnom(Kep2(6), Kep2(2)) - KepChange(5);
KepF(6)     = E - Kep2(2)*sin(E);
Y_F         = State(const.GM_Earth,KepF,0);


Y_F = State(const.GM_Earth,KepF,0);
hold on
plot3(Y_F(1),Y_F(2),Y_F(3),'x')
hold on
plot3(Y_S(length(Y_S),1),Y_S(length(Y_S),2),Y_S(length(Y_S),3),'x')

Y2  = State(const.GM_Earth,Kep(2,:),-((Mjd_Epoch(2)-Mjd_Epoch(1))*86400 - Eph(length(Eph),1)));
hold on
plot3(Y2(1),Y2(2),Y2(3),'or')
Error = ((Y_F' - Y2').^2)';
sqrt(sum(Error(1:3,:)))'
ans/((Mjd_Epoch(2)-Mjd_Epoch(1))*86400)
sqrt(sum(Error(4:6,:)))'






