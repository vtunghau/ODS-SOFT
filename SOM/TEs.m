%% Update State
Y_S         = State(const.GM_Earth,Kep(1,:),Eph(length(Eph),1));
Kep3        = Element(const.GM_Earth,Eph(length(Eph),2:7))
Kep2        = Element(const.GM_Earth,Y_S);
KepChange   = Kep3 - Kep2;
KepF        = Kep3;

E           = EccAnom(Kep2(6), Kep2(2));
v           = 2*atan(sqrt((1+Kep2(2))/(1-Kep2(2)))*tan(E/2));
v           = 2*pi + v - KepChange(5) ;
KepF(6)     = mod(v2M(v,Kep2(2)),2*pi) + (TD.TD1(1)/86400)*Sample^2;
KepF(6)     = mod(v2M(v,Kep2(2)),2*pi) - sin(Kep2(3))*KepChange(4);
Y_F         = State(const.GM_Earth,KepF,0);



Y2  = State(const.GM_Earth,Kep(2,:),-((Mjd_Epoch(2)-Mjd_Epoch(1))*86400 - Eph(length(Eph),1)));
Y2  = State(const.GM_Earth,Kep(2,:),-0.018741131);
Error           = ((Y_F' - Y2').^2)';
Error_dis       = sqrt(sum(Error(1:3,:)))'
Error_dis_dt    = Error_dis/((Mjd_Epoch(2)-Mjd_Epoch(1))*86400);
Error_vel       = sqrt(sum(Error(4:6,:)))';

plot3(Y_F(1),Y_F(2),Y_F(3),'or');               hold on
plot3(Y2(1),Y2(2),Y2(3),'og');                  hold on

