%% Making a photo contain a satellite path in the sky
% input
% P1: point 1
% P2: point 2
% pixel size


function out = path2im(P1, P2, pixelhz,pixelver)
S1  = [P1(1); P2(1)];
S2  = [P1(2); P2(2)];
COE = [S1(1) 1;S1(2) 1]\S2;

FWHMY = 5/pixelhz;
FWHMX = sqrt((S1(1)-S1(2))^2+(S2(1)-S2(2))^2)/2;
k1    = 50;
k2    = 1;
k3    = 3;

wx  = FWHMX/(log(2)^(1/(k1*k3)));
wy  = FWHMY/(log(2)^(1/(k2*k3)));
out   = zeros(pixelver,pixelhz);
for i = 1:pixelver
    for j = 1:pixelhz
        xi =  (i/(pixelver/2) - 1)+(S2(1)+S2(2))/2;
        yi =  (j/(pixelhz/2) - 1)-(S1(1)+S1(2))/2;
        x  = xi.*cos(-(pi/2 + atan(COE(1)))) - yi.*sin(-(pi/2 + atan(COE(1))));
        y  = xi.*sin(-(pi/2 + atan(COE(1)))) + yi.*cos(-(pi/2 + atan(COE(1))));
        out(i,j) = exp(-(abs(x/wx).^k1+abs(y/wy).^k2).^k3);
    end
end
