addpath('..\SGP4\');


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
                   
[satrec, ro ,vo] = sgp4 (satrec,  0.0);

tsince = startmfe;

jd = satrec.jdsatepoch + tsince/1440.0
    
[year,mon,day,hr,minute,sec] = invjday ( jd );

app.EpochTime.Value                = sprintf('%02d/%02d/%04d - %02d:%02d:%02d',day,mon,year,hr,minute,round(sec));
%               app.EpochTime              
% [p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);

% app.year.Value  = year;
% app.month.Value = mon;
% app.day.Value   = day;
% app.hour.Value  = hr;
% app.min.Value   = minute;
% app.sec.Value   = sec;


% lon_Sta         = app.lon_Sta.Value*const.Rad;          % [rad]
% lat_Sta         = app.lat_Sta.Value*const.Rad;          % [rad]
% alt_h           = app.alt_h.Value;                      % [m]

lon_Sta         = 0;          % [rad]
lat_Sta         = 0;          % [rad]
alt_h           = 0;                      % [m]

R_Sta           = Position(lon_Sta, lat_Sta, alt_h)/1000;    % Geocentric position vector
E               = LTCMatrix(lon_Sta, lat_Sta);          % Transformation to local tangent coordinates

        
%% Checking Visible Time on the sky
mm = 1;
mn = 0;
for i = 1:5000
    
    tsince = tsince + 1;
    
    [satrec, ro, vo] = sgp4 (satrec,  tsince);
    
    jd = satrec.jdsatepoch + tsince/1440.0;
    
    [year,mon,day,hr,minute,sec] = invjday ( jd );
    
    Mjd_UTC                                     = Mjday(year,mon,day,hr,minute,sec);
    
    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
    
    Mjd_UT1                                     = Mjd_UTC + UT1_UTC/86400;
    
    U           = R_z(gmst(Mjd_UT1));                   % Earth rotation 
            
    s           = E * ( U*ro' - R_Sta' );               % Topocentric position vector

    Dist = norm(s);                                     % Distance
    
    [Azim, Elev]            = AzEl(s);                  % Azimuth, Elevation
    
%     Azel_Data(i,1)    = Mjd_UTC;
%     
%     Azel_Data(i,2:3)  = [Azim, Elev];
%     
%     Azel_Data(i,4:9)  = [ro vo];
%     
%     Azel_Data(i,10)   = Dist;
    
    if (Elev > 0)&(mn == 0)
        Vis(mm,1)     = jd;
        fprintf('%.10f \n',Elev);
        mn              = 1;
    end
    
    if (Elev < 0)&(mn == 1)
        Vis(mm,2)    = jd;
        Vis(mm,3)      = (Vis(mm,2) - Vis(mm,1))*1440;
        app.VisibleList.Items{end + 1} = sprintf('Time %03d',mm);
        mm             = mm + 1;
        mn              = 0;
    end
    
end

save ..\Data\Vis.mat Vis

rmpath('..\SGP4\');









