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
AuxParam.n          = 70;
AuxParam.m          = 70;
AuxParam.sun        = 1;
AuxParam.moon       = 1;
AuxParam.planets    = 1;
AuxParam.sRad       = 0;
AuxParam.drag       = 0;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;


ge = 398600.4418; % Earth gravitational constant
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

% TLE file name
fname = 'tle.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

% 19-32	04236.56031392	Element Set Epoch (UTC)
% 3-7	25544	Satellite Catalog Number
% 9-16	51.6335	Orbit Inclination (degrees)
% 18-25	344.7760	Right Ascension of Ascending Node (degrees)
% 27-33	0007976	Eccentricity (decimal point assumed)
% 35-42	126.2523	Argument of Perigee (degrees)
% 44-51	325.9359	Mean Anomaly (degrees)
% 53-63	15.70406856	Mean Motion (revolutions/day)
% 64-68	32890	Revolution Number at Epoch

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
a = ( ge/(no*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
rNo = str2num(tline(65:end));               % Revolution Number at Epoch

fclose(fid);

%% SGP4 initial propertise
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
satdata.xndt2o = 0 * TD1 * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o = 0 * TD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar = BStar*0;
satdata_update = satdata;

%% Orbital Probagation Methodology

Y0   = Y(1,:)';
time = 0;
KEP  = Kep(1,:);
Up   = 1000;
TIME = (Mjd_Epoch(2)-Mjd_Epoch(1))*86400;
N    = fix(TIME/Up)+1;
for counter = 1:N
    time = Up*counter
    if(time > TIME)
        Sample  = fix(mod(TIME-time-Up,Up));
    else
        Sample  = Up;
    end
    %% Update Orbit Element for Ephemeris
    Eph         = Ephemeris(Y0,Sample,Step);
    lengthEph   = size(Eph);
    lengthEph   = lengthEph(1);
    [rte, vte]  = sgp4((time - Up + Sample)/60, satdata);
    Y_S         = [rte; vte]*1000;
    Kep3        = Element(const.GM_Earth,Eph(lengthEph,2:7));
    Kep2        = Element(const.GM_Earth,Y_S);
    E           = EccAnom(Kep2(6), Kep2(2));
    v           = 2*atan(sqrt((1+Kep2(2))/(1-Kep2(2)))*tan(E/2));
    KepChange   = (Kep3 - Kep2)';
    KEP         = Kep3;
    KEP(1)      = Kep2(1);
    if(Kep3(5) > Kep2(5))
        v       = 2*pi + v - KepChange(5);
        KEP(6)  = mod(v2M(v,Kep2(2)),2*pi);
    else
        v       = v - KepChange(5);
        KEP(6)  = mod(v2M(v,Kep2(2)),2*pi);
    end

    %% Correct State
    Y0          = (State(const.GM_Earth,KEP,0));
    ErrorR      = Y0 - Eph(lengthEph,2:7)'

    %% Update Element
    AuxParam.Mjd_UTC       = AuxParam.Mjd_UTC + Eph(lengthEph,1)/86400;    
    plot3(Eph(:,2),Eph(:,3),Eph(:,4),'k');hold on

end


%% Check accuracy of prediction
Y0              = (State(const.GM_Earth,KEP,((Mjd_Epoch(2)-Mjd_Epoch(1))*86400 - (time - Up + Sample))))';
Y2              = Y(2,:);
Error           = ((Y0 - Y2).^2)';
Error_dis       = sqrt(sum(Error(1:3,:)))'
Error_dis_dt    = Error_dis/((Mjd_Epoch(2)-Mjd_Epoch(1))*86400);
Error_vel       = sqrt(sum(Error(4:6,:)))';
plot3(Y0(1),Y0(2),Y0(3),'xr');                  hold on
plot3(Y2(1),Y2(2),Y2(3),'or');                  hold on
% %zoomPoint(Y2,100000)
[Xx,Yx,Zx] = sphere(50);
Xx = Xx*const.R_Earth;
Yx = Yx*const.R_Earth;
Zx = Zx*const.R_Earth;
surf(Xx,Yx,Zx)
[rte, vte]      = sgp4(((Mjd_Epoch(2)-Mjd_Epoch(1))*86400)/60, satdata);
Y0              = [rte; vte]*1000;
Error           = ((Y0' - Y2).^2)';
Error_dis       = sqrt(sum(Error(1:3,:)))'
