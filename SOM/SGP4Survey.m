%% Satellite probagation
TLE2RV
clc
%close all
clear all
%close all
format long g
global AuxParam eopdata PC Cnm Snm const swdata

%% Load initial data

SAT_Const
load DE430Coeff.mat
load Epoch.mat
load State.mat
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


ge = 398600.4418; % Earth gravitational constant
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

% TLE file name
fname = 'tle.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

% read first line
tline = fgetl(fid);
Cnum = tline(3:7);      			        % Catalog Number (NORAD)
SC   = tline(8);					        % Security Classification
ID   = tline(10:17);			            % Identification Number
year = str2num(tline(19:20));               % Year
doy  = str2num(tline(21:32));               % Day of year
epoch = str2num(tline(19:32));              % Epoch
TD1   = str2num(tline(34:43));              % first time derivative
TD2   = str2num(tline(45:50));              % 2nd Time Derivative
ExTD2 = tline(51:52);                       % Exponent of 2nd Time Derivative
BStar = str2num(tline(54:59));              % Bstar/drag Term
ExBStar = str2num(tline(60:61));            % Exponent of Bstar/drag Term
BStar = BStar*1e-5*10^ExBStar;
Etype = tline(63);                          % Ephemeris Type
Enum  = str2num(tline(65:end));             % Element Number

% read second line
tline = fgetl(fid);
i = str2num(tline(9:16));                   % Orbit Inclination (degrees)
raan = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
e = str2num(strcat('0.',tline(27:33)));     % Eccentricity
omega = str2num(tline(35:42));              % Argument of Perigee (degrees)
M = str2num(tline(44:51));                  % Mean Anomaly (degrees)
no = str2num(tline(53:63));                 % Mean Motion
a = ( ge/(no*2*pi/86400)^2 )^(1/3)         % semi major axis (m)
rNo = str2num(tline(65:end));               % Revolution Number at Epoch

fclose(fid);

satdata.epoch = epoch;
satdata.norad_number = Cnum;
satdata.bulletin_number = ID;
satdata.classification = SC; % almost always 'U'
satdata.revolution_number = rNo;
satdata.ephemeris_type = Etype;
satdata.xmo = M * (pi/180);
satdata.xnodeo = raan * (pi/180);
satdata.omegao = omega * (pi/180);
satdata.xincl = i * (pi/180);
satdata.eo = e;
satdata.xno = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o = 0*TD1 * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o = 0*TD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar = 0*BStar;

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

[rte, vte]  = sgp4(0, satdata);
Y_S         = [rte; vte]*1000;
EPH         = Ephemeris(Y_S,((Mjd_Epoch(2) - Mjd_Epoch(1))*1440+1),60);

for i = 1:((Mjd_Epoch(2) - Mjd_Epoch(1))*1440+2)
    [rte, vte]  = sgp4(i-1, satdata);
    Y_S(i,1:6)         = ([rte; vte]*1000)';
%     KEP(i,:)    = Element(const.GM_Earth,Y_S);
end

% for i = 1:length(EPH)
%     KE(i,:)  = Element(const.GM_Earth,EPH(i,2:7));
% end
% 
% close all
% figure
% grid on
% plot(KEP(:,6)+KEP(:,5))
% hold on
% % plot(KE(:,6)+KE(:,5))
% for i = 1:10
%     plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,6)+Kep(i,5),'xr')  
% end
% %zoomDATA([(Mjd_Epoch(2) - Mjd_Epoch(1))*1440+1 Kep(2,6)+Kep(2,5)],0.1);
% figure
% plot(KEP(:,1))
% hold on
% for i = 1:10
%     plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,1),'xr')  
% end
% ylim([(a)*1000-10000 (a)*1000+10000])


