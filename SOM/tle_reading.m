app.Lamp.Color = 'blue';
app.Lamp.Enable = 'on'
try
    format long g
    global eopdata const

    %% Load initial data
    SAT_Const
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
    fname = app.Tle_dir.Value;
    fid = fopen(fname, 'r');

    % read first line
    tline = fgetl(fid);
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID   = tline(10:17);			            % Identification Number
    year = str2num(tline(19:20)) ;              % Year
    doy  = str2num(tline(21:32)) ;              % Day of year
    epoch = str2num(tline(19:32)) ;             % Epoch
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
    i       = str2num(tline(9:16));                % Orbit Inclination (degrees)
    raan    = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e       = str2num(strcat('0.',tline(27:33)));  % Eccentricity
    omega   = str2num(tline(35:42));               % Argument of Perigee (degrees)
    M       = str2num(tline(44:51));               % Mean Anomaly (degrees)
    no      = str2num(tline(53:63));               % Mean Motion
    a       = ( ge/(no*2*pi/86400)^2 )^(1/3);      % semi major axis (m)
    rNo     = str2num(tline(65:end));              % Revolution Number at Epoch

    fclose(fid);

    satdata.epoch = epoch;
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
    if (year < 57)
        year = year + 2000;
    else
        year = year + 1900;
    end
    [mon,day,hr,minute,sec] = days2mdh(year,doy);
    Mjd_Epoch = Mjday(year,mon,day,hr,minute,sec);

    [rte, vte]  = sgp4(0, satdata);
    Y           = [rte; vte]*1000;
    Kep         = Element(const.GM_Earth,Y);

    app.semi.Value = Kep(1);
    app.ecc.Value  = Kep(2);
    app.inc.Value  = Kep(3);
    app.raad.Value = Kep(4);
    app.omega.Value= Kep(5);
    app.MeanN.Value= Kep(6);

    app.year.Value  = year;
    app.month.Value = mon;
    app.day.Value   = day;
    app.hour.Value  = hr;
    app.min.Value   = minute;
    app.sec.Value   = sec;
    app.Mjd_epoch.Value = Mjd_Epoch;
    app.Status.Text = 'Reading TLE successfully'; 
    app.Lamp.Color = 'green';
    app.Log.Value  = '';
    
catch ME
    app.Lamp.Color = 'red';
    ME
    app.Log.Value  = ME.message;
end
