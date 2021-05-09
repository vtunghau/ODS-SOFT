%% Satellite probagation
TLE2RV
clc
clear all
format long g
global eopdata const

%% Load initial data

SAT_Const

%% Load File 
load ..\Data\Input.mat
% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);


ge = 398600.4418; % Earth gravitational constant
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

% TLE file name
fname = Input.tle_path;

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
satdata.xndt2o = TD1 * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o = TD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar = BStar;

%% Setup for Az El Probagation

lon_Sta     = Input.lon_Sta*const.Rad; % [rad]
lat_Sta     = Input.lat_Sta*const.Rad; % [rad]
alt_h       = Input.alt_h;     % [m]
Step_size   = Input.Step; % second
N_obs       = Input.N_Step;
R_Sta       = Position(lon_Sta, lat_Sta, alt_h); % Geocentric position vector
E           = LTCMatrix(lon_Sta, lat_Sta);       % Transformation to local tangent coordinates

if (year < 57)
    year = year + 2000;
else
    year = year + 1900;
end
[mon,day,hr,minute,sec] = days2mdh(year,doy);
Mjd_Epoch = Mjday(year,mon,day,hr,minute,sec);
Mjd_Epoch = fix(Mjd_Epoch) + mod(epoch,1);
counter = 1;
inn     = 1;
while counter <= N_obs
    [rte, vte]  = sgp4((inn-1)*Step_size/60, satdata);
    Y         = [rte; vte]*1000;
    
    Mjd_UTC     = Mjd_Epoch + (inn-1)*Step_size/86400;
    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
    Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
    dt = (Mjd_UTC-Mjd_Epoch)*86400;        % Time since epoch [s] 

    r = Y(1:3);
    U = R_z(gmst(Mjd_UT1));                % Earth rotation 
    s = E * ( U*r - R_Sta' );              % Topocentric position vector

    Dist = norm(s);                        % Distance
    [Azim, Elev]        = AzEl(s);                % Azimuth, Elevation
    if(Elev > )
        (inn-1)*Step_size
        Eph(counter,:)  = Y';
        Azel_Data(counter,1)    = Mjd_UTC;
        Azel_Data(counter,2:3)  = [Azim, Elev];
        counter = counter + 1;
    end
    inn = inn + 1;
end
% for inn = 1:N_obs
%     [rte, vte]  = sgp4((inn-1)*Step_size/60, satdata);
%     Y         = [rte; vte]*1000;
%     Eph(inn,:)  = Y';
%     Mjd_UTC     = Mjd_Epoch + (inn-1)*Step_size/86400;
%     [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
%     Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
%     dt = (Mjd_UTC-Mjd_Epoch)*86400;        % Time since epoch [s] 
% 
%     r = Y(1:3);
%     U = R_z(gmst(Mjd_UT1));                % Earth rotation 
%     s = E * ( U*r - R_Sta' );              % Topocentric position vector
% 
%     Dist = norm(s);                        % Distance
%     [Azim, Elev]        = AzEl(s);                % Azimuth, Elevation
%     Azel_Data(inn,1)    = Mjd_UTC;
%     Azel_Data(inn,2:3)  = [Azim, Elev];
% end

Inertial.R_Sta = R_Sta;
Inertial.E     = E;
save ..\Data\Azel_Data.mat Azel_Data
save ..\Data\Ephemeris.mat Eph
save ..\Data\Inertial.mat Inertial