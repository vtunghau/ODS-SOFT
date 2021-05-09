satdata.epoch = app.Mjd_epoch.Value;
%satdata.bulletin_number = ID;
%satdata.classification = SC; % almost always 'U'
%satdata.revolution_number = rNo;
%satdata.ephemeris_type = Etype;
ge = 398600.4418; % Earth gravitational constant
satdata.xmo = app.MeanN.Value;
satdata.xnodeo = app.raad.Value;
satdata.omegao = app.omega.Value;
satdata.xincl = app.inc.Value;
satdata.eo = app.ecc.Value;
satdata.xno = sqrt(ge/app.semi.Value^3)/(2*pi/86400);
satdata.xndt2o = 0;
satdata.xndd6o = 0;
satdata.bstar = 0;
Step_size     = 10;

for inn = 1:50
    [rte, vte]  = sgp4((inn-1)*Step_size/60, satdata);
    Y         = [rte; vte]*1000;
    Eph(inn,:)  = Y';
end
app.Status.Text = 'Propagating the orbital object successfully'; 
