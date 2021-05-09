app.Lamp.Enable = 'on';
app.Lamp.Color  = 'blue';
app.Status.Text = 'Propagating the orbital object ...'; 
app.azelList.Items = {};
try 
    addpath('..\SGP4\');

    load ..\Data\Vis.mat
    % initial setup

    % global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  

    global opsmode eopdata const 

    fid = fopen('eop19620101.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
    fclose(fid);

    SAT_Const

    opsmode     = 'a';

    % -- init --
    whichconst  = 721;
    rad         = 180.0 / pi;
    infilename  = app.Tle_dir.Value;
    typerun     = 'c';
    infile      = fopen(infilename, 'r');
    longstr1    = fgets(infile, 130);
    longstr2    = fgets(infile, 130);
    typeinput   = 'e';
    
    [satrec, startmfe, stopmfe, deltamin] = twoline2rv( whichconst, ...
                           longstr1, longstr2, typerun, typeinput);
    
%    [satrec, ro ,vo]  = sgp4(satrec,  );
%    satrec.jdsatepoch = Vis(1,1);
    
    
    value           = app.VisibleList.Value;
    itemselection   = char(value);
    tsince      = (Vis(str2num(itemselection(6:8)),1) - satrec.jdsatepoch)*86400;

    lon_Sta         = app.lon_Sta.Value*const.Rad;          % [rad]
    lat_Sta         = app.lat_Sta.Value*const.Rad;          % [rad]
    alt_h           = app.alt_h.Value;                      % [m]
    
    Sta.lon       = lon_Sta;
    Sta.lat       = lat_Sta;
    Sta.alt       = alt_h;

%     lon_Sta         = 0;          % [rad]
%     lat_Sta         = 0;          % [rad]
%     alt_h           = 0;                      % [m]

    R_Sta           = Position(lon_Sta, lat_Sta, alt_h)/1000;    % Geocentric position vector
    E               = LTCMatrix(lon_Sta, lat_Sta);          % Transformation to local tangent coordinates
    Inertial.R_Sta  = R_Sta;
    Inertial.E      = E;
    save ..\Data\Inertial.mat Inertial
    fprintf('%.10f \n',satrec.jdsatepoch);
    
    for i = 1:app.noo.Value
        app.azelList.Items{end + 1} = sprintf('Item %03d',i);
        % Sec 1

        [satrec, ro, vo] = sgp4 (satrec,  tsince/60);

        jd = satrec.jdsatepoch + tsince/86400.0;

        [year,mon,day,hr,minute,sec] = invjday ( jd );
        
        Prop_list(i,1)               = jd;

        Mjd_UTC                                     = Mjday(year,mon,day,hr,minute,sec);

        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');

        Mjd_UT1                                     = Mjd_UTC + UT1_UTC/86400;

        U           = R_z(gmst(Mjd_UT1));                   % Earth rotation 

        s           = E * ( U*ro' - R_Sta' );               % Topocentric position vector

        Dist = norm(s);                                     % Distance

        [Azim, Elev]            = AzEl(s);                  % Azimuth, Elevation

        Azel_Data(i*3-2,1)    = Mjd_UTC;
        Azel_Data(i*3-2,2:3)  = [Azim, Elev];
        Azel_Data(i*3-2,4:9)  = [ro vo];
        Azel_Data(i*3-2,10)   = Dist;
        [Az, El]              = AzEl2RaDec(Azim*const.Deg,Elev*const.Deg,lat_Sta,lon_Sta, Mjd_UT1);
        Azel_Data(i*3-2,11:12)= [Az*const.Rad, El*const.Rad] ;
        % Sec 2

        tsince = tsince + app.obs_time.Value/2;

        [satrec, ro, vo] = sgp4 (satrec,  tsince/60);

        jd = satrec.jdsatepoch + tsince/86400.0;

        [year,mon,day,hr,minute,sec] = invjday ( jd );

        Mjd_UTC                                     = Mjday(year,mon,day,hr,minute,sec);

        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');

        Mjd_UT1                                     = Mjd_UTC + UT1_UTC/86400;

        U           = R_z(gmst(Mjd_UT1));                   % Earth rotation 

        s           = E * ( U*ro' - R_Sta' );               % Topocentric position vector

        Dist = norm(s);                                     % Distance

        [Azim, Elev]            = AzEl(s);                  % Azimuth, Elevation
        
        Prop_list(i,3:4)               = [Azim, Elev] ;

        Azel_Data(i*3-1,1)    = Mjd_UTC;
        Azel_Data(i*3-1,2:3)  = [Azim, Elev];
        Azel_Data(i*3-1,4:9)  = [ro vo];
        Azel_Data(i*3-1,10)   = Dist;
        
        [Az, El]                = AzEl2RaDec(Azim*const.Deg,Elev*const.Deg,lat_Sta,lon_Sta, Mjd_UT1);
        Azel_Data(i*3-1,11:12)  = [Az*const.Rad, El*const.Rad] ;
        Azel_Data(i*3-1,13)     = Mjd_UT1;
        AZ(i,:) = [Mjd_UTC,0,0,Azim,Elev,ro,vo,Dist];
        fprintf('%f \n',Mjd_UTC);
        % Sec 3

        tsince = tsince + app.obs_time.Value/2;

        [satrec, ro, vo] = sgp4 (satrec,  tsince/60);

        jd = satrec.jdsatepoch + tsince/86400.0;
        
        Prop_list(i,2)               = jd;

        [year,mon,day,hr,minute,sec] = invjday ( jd );

        Mjd_UTC                                     = Mjday(year,mon,day,hr,minute,sec);

        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');

        Mjd_UT1                                     = Mjd_UTC + UT1_UTC/86400;

        U           = R_z(gmst(Mjd_UT1));                   % Earth rotation 

        s           = E * ( U*ro' - R_Sta' );               % Topocentric position vector

        Dist = norm(s);                                     % Distance

        [Azim, Elev]            = AzEl(s);                  % Azimuth, Elevation

        tsince              = tsince - app.obs_time.Value/2 + app.step_time.Value/2;
        Azel_Data(i*3,1)    = Mjd_UTC;
        Azel_Data(i*3,2:3)  = [Azim, Elev];
        Azel_Data(i*3,4:9)  = [ro vo];
        Azel_Data(i*3,10)   = Dist;
        [Az, El]            = AzEl2RaDec(Azim*const.Deg,Elev*const.Deg,lat_Sta,lon_Sta, Mjd_UT1);
        Azel_Data(i*3,11:12)= [Az*const.Rad, El*const.Rad] ;

    end
    save ..\Data\Prop_list.mat Prop_list
    save ..\Data\Azel_Data.mat Azel_Data
    save ..\Data\AZ.mat AZ
    save ..\Data\Sta.mat Sta
    app.Status.Text     = 'Propagate the orbital object successfully'; 
        
    app.Lamp.Color      = 'green';
        
    app.Log.Value       = '';
        
catch ME
    app.Status.Text     = 'Propagate the orbital object failed'; 
        
    app.Lamp.Color      = 'red';
        
    app.Log.Value       = ME.message;
end

rmpath('..\SGP4\');