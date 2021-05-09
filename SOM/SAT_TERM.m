close all
clear all

%% SGP4 initial propertise
satdata.epoch = epoch;
satdata.norad_number = 'NULL';
satdata.bulletin_number = 'NULL';
satdata.classification = 'U'; % almost always 'U'
satdata.revolution_number = 0;
satdata.ephemeris_type = '0';
satdata.xmo = M * (pi/180);
satdata.xnodeo = raan * (pi/180);
satdata.omegao = omega * (pi/180);
satdata.xincl = i * (pi/180);
satdata.eo = e;
satdata.xno = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o = 0;
satdata.xndd6o = 0 ;
satdata.bstar = 0;
