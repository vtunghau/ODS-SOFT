clc
close all
nuy1    = sqrt(const.GM_Earth/(Kep(1,1)^3))
T1      = 2*pi/nuy1 
nuy2    = sqrt(const.GM_Earth/(Kep(2,1)^3))
T2      = 2*pi/nuy2
Dnuy    = (nuy2 - nuy1)/55048;
t       = 1:length(Eph);
Nuy     = Dnuy.*t;
for i = 1:length(Eph)
    KE(i,:)  = Element(const.GM_Earth,Eph(i,2:7));
    nuy(i,:) = sqrt(const.GM_Earth/KE(i,1)^3);
    T(i,:)   = 2*pi/nuy(i,:);
end

figure 
plot(nuy); hold on
plot(length(Eph),nuy2,'x'); hold on


DT   = (T' + (0.005/86400.*t));
nnuy = 2*pi./DT;
plot(nnuy);
Ap   = (const.GM_Earth./((nnuy).^2)).^(1/3);

figure 
plot(Ap);hold on
plot(55048,Kep(2,1),'x'); hold on
plot(KE(:,1))