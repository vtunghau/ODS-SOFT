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
infilename  = 'SGP4-VER.TLE';
typerun     = 'c';
infile      = fopen(infilename, 'r');

for i = 1:50
    longstr1(i,:)    = fgets(infile, 130);
    longstr2(i,:)    = fgets(infile, 130);
end
typeinput   = 'e';

[satrecP, ~, ~, ~] = twoline2rv( whichconst, ...
                       longstr1(1,:), longstr2(1,:), typerun, typeinput);
satreOLD =           satrecP;   

for i = 1:length(longstr1(:,1))-1
    [satrec, ~, ~, ~] = twoline2rv( whichconst,longstr1(i+1,:), longstr2(i+1,:), typerun, typeinput);
    
    if satreOLD.jdsatepoch < satrec.jdsatepoch
                   
        fprintf('%.5f \n',satrec.jdsatepoch - satrecP.jdsatepoch);

        [~, ro, ~]     = sgp4 (satrec,  0);

        [~, roP, ~]    = sgp4 (satrecP,  (satrec.jdsatepoch - satrecP.jdsatepoch)* 1440);

        satrecD(i,:)        = [satrec.jdsatepoch sqrt(sum((roP - ro).^2))];  

        satreOLD       = satrec;
    end
end
                   

rmpath('..\SGP4\');









