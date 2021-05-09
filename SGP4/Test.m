clc;
close all;
clear all;

global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  

global opsmode

opsmode     = 'a';

% -- init --
whichconst  = 721;
rad         = 180.0 / pi;
infilename  = 'SGP4-VER.TLE';
typerun     = 'c';
infile      = fopen(infilename, 'r');
longstr1    = fgets(infile, 130);
longstr2    = fgets(infile, 130);
typeinput   = 'e';

[satrec, startmfe, stopmfe, deltamin] = twoline2rv( whichconst, ...
                       longstr1, longstr2, typerun, typeinput);
                   
[satrec, ro ,vo] = sgp4 (satrec,  0.0);

tsince = startmfe;

jd = satrec.jdsatepoch + tsince/1440.0;

[year,mon,day,hr,minute,sec] = invjday ( jd );

[p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);

% if ( abs(tsince) > 1.0e-8 )
%                 tsince = tsince - deltamin;
% end
% 
% i = 1;
% while ((tsince < stopmfe) && (satrec.error == 0))
% 
%     tsince = tsince + deltamin
% 
%     if(tsince > stopmfe)
% 
%         tsince = stopmfe;
% 
%     end
%     
%     [satrec, ro, vo] = sgp4 (satrec,  tsince);
%     
%     jd = satrec.jdsatepoch + tsince/1440.0;
%     
%     [year,mon,day,hr,minute,sec] = invjday ( jd );
%                                  
%     [p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);
%     
%     result(i,:) = [year,mon,day,hr,minute,sec,p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper, ro, vo];
%     
%     i = i + 1;
% end