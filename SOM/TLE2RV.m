clc
clear
format long g

SAT_Const

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

fid = fopen('TLE.txt', 'r');
for counter = 1:10
    tline = fgetl(fid);
    Cnum(counter,:) = tline(3:7);      			        % Catalog Number (NORAD)
    SC(counter,:)   = tline(8);					        % Security Classification
    ID(counter,:)   = tline(10:17);			            % Identification Number
    year(counter) = str2num(tline(19:20));               % Year
    doy(counter)  = str2num(tline(21:32));               % Day of year
    epoch(counter) = str2num(tline(19:32));              % Epoch
    TD1(counter)   = str2num(tline(34:43));              % first time derivative
    TD2(counter)   = str2num(tline(45:50));              % 2nd Time Derivative
    ExTD2(counter,:) = tline(51:52);                       % Exponent of 2nd Time Derivative
    BStar(counter) = str2num(tline(54:59));              % Bstar/drag Term
    ExBStar(counter) = str2num(tline(60:61));            % Exponent of Bstar/drag Term
    BStar(counter) = BStar(counter)*1e-5*10^ExBStar(counter);
    Etype(counter,:) = tline(63);                          % Ephemeris Type
    Enum(counter)  = str2num(tline(65:end));             % Element Number

    % read second line
    tline = fgetl(fid);
    i(counter) = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    raan(counter) = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e(counter) = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    omega(counter) = str2num(tline(35:42));              % Argument of Perigee (degrees)
    M(counter) = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no(counter) = str2num(tline(53:63));                 % Mean Motion
    a(counter) = ( const.GM_Earth/10e8/(no(counter)*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo(counter) = str2num(tline(65:end));               % Revolution Number at Epoch
end
%load Kep.mat
%load Epoch.mat
% load GPS.mat
% load GPSepoch.mat
% Kep       = Orbit([1,2],1:6);
% Mjd_Epoch = Epoch([1,2]);
for counter = 1:10
    satdata.epoch = epoch(counter);
    satdata.norad_number = Cnum(counter,:);
    satdata.bulletin_number = ID(counter,:);
    satdata.classification = SC(counter,:); % almost always 'U'
    satdata.revolution_number = rNo(counter);
    satdata.ephemeris_type = Etype(counter,:);
    satdata.xmo = M(counter) * (pi/180);
    satdata.xnodeo = raan(counter) * (pi/180);
    satdata.omegao = omega(counter) * (pi/180);
    satdata.xincl = i(counter) * (pi/180);
    satdata.eo = e(counter);
    satdata.xno = no(counter) * TWOPI / MINUTES_PER_DAY;
    satdata.xndt2o = TD1(counter) * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
    satdata.xndd6o = TD2(counter) * TWOPI / MINUTES_PER_DAY_CUBED;
    satdata.bstar = BStar(counter);
    year(counter) = year(counter) + 2000;
    [rte, vte]  = sgp4(0, satdata);
    Y(counter,:)         = ([rte; vte]*1000)';
    Kep(counter,:)    = Element(const.GM_Earth,Y(counter,:)'); 
    
    [mon,day,hr,minute,sec] = days2mdh(year(counter),doy(counter));
     
    Mjd_Epoch(counter) = Mjday(year(counter),mon,day,hr,minute,sec);

    %Kep(i,:) = [a(i),e(i),inc(i),Omega(i),omega(i),M0(i)]; % Keplerian elements

    %Y(i,:) = State(const.GM_Earth, Kep(i,:), 0);          % Inertial vector
    %[r,v] = orb2rv(p(i)/1000,e(i),inc(i),Omega(i),omega(i),E(i));
    %Y(i,:) = [r;v]'*1000;
end

save Epoch.mat Mjd_Epoch
save State.mat Y
save Kep.mat Kep
