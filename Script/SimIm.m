clc;
close all;
clear all;

% Featured for Database
% Input
fid     = fopen('sao','r');
YO      = 2020; 
i       = 0;
FOV     = 5;
minStar = 8;
minPair = 8;
azi     = 5;
ele     = 40;
ori     = 45;
azi     = azi*15;
intensity = 1;


%setup
YYYY    = YO - 2000;
SAO     = zeros(258697,7);




while (~feof(fid))
    tline = fgetl(fid);
    if(str2double(tline(81:84))<=minStar)
        i = i+1;
        SAO(i,1)= str2double(tline(1:6));    % ID
        SAO(i,2)= str2double(tline(81:84));  %   81- 84  F4.1   mag     Vmag     []?=99.9 Visual magnitude
        COLOR   = tline(85);                 %  Color
        SAO4    = str2double(tline(151:152));%  151-152  I2     h       RA2000h  Hours RA, equinox, epoch J2000.0
        SAO5    = str2double(tline(153:154));%  153-154  I2     min     RA2000m  Minutes RA, equinox, epoch J2000.0
        SAO6    = str2double(tline(155:160));%  155-160  F6.3   s       RA2000s  Seconds RA, equinox, epoch J2000.0
        SAO7    = str2double(tline(161:167));%  161-167  F7.4   s/a     pmRA2000 Annual proper motion in FK5 system
        SAO8    = double(tline(168)==43)*2-1;             %      168  A1     ---     DE2000-  Sign Dec, equinox, epoch J2000.0
        SAO9    = str2double(tline(169:170));%  169-170  I2     deg     DE2000d  Degrees Dec, equinox, epoch J2000.0
        SAO10   = str2double(tline(171:172));%  171-172  I2     arcmin  DE2000m  Minutes Dec, equinox, epoch J2000.0
        SAO11   = str2double(tline(173:177));%  173-177  F5.2   arcsec  DE2000s  Seconds Dec, equinox, epoch J2000.0
        SAO12   = str2double(tline(178:183));%  178-183  F6.3  arcsec/a pmDE2000 ? Annual proper motion in FK5 system (10)

        SAO(i,3)= (SAO4 + SAO5/60 + SAO6/3600 + SAO7*YYYY/3600)*15;
        SAO(i,4)= (SAO9 + SAO10/60 + SAO11/3600 + SAO12*YYYY/3600)*SAO8;
        
        switch COLOR
            case 79
                SAO(i,5:7) = [146 181 255]/255;
            case 66
                SAO(i,5:7) = [170 191 255]/255;
            case 65
                SAO(i,5:7) = [202 212 255]/255;
            case 70
                SAO(i,5:7) = [248 247 255]/255;
            case 71
                SAO(i,5:7) = [255 244 234]/255;
            case 75
                SAO(i,5:7) = [255 210 161]/255;
            case 77
                SAO(i,5:7) = [255 204 111]/255;
            otherwise
                SAO(i,5:7) = [255 255 255]/255;
        end
    end
end
SAO = SAO(SAO(:,1)>0,:);
save SAODATA.mat SAO;
SAO = SAO(SAO(:,2)<minPair,:);

hFOV = FOV/2;
SAO(:,8)    =  2*asind(sqrt(sind((SAO(:,4)-ele)/2).^2 + cosd(SAO(:,4)).*cosd(ele).*(sind((SAO(:,3)-azi)/2).^2)));
SAO         = SAO(SAO(:,8)<hFOV,:);



SAO         = SAO(SAO(:,1)>0,:);
[~, index]  = sort(SAO(:,2),'ascend');
SAO         = SAO(index,1:7);
length      = length(SAO);
Pair        = zeros(10e6,6);
m           = 0;

for i = 1:length-2
    for j = i+1:length-1
        p1 = angular_distance(SAO(i,4:-1:3),SAO(j,4:-1:3));
        if  p1 < FOV
            for k = j+1:length
                p2 = angular_distance(SAO(i,4:-1:3),SAO(k,4:-1:3));
                if p2 <  FOV
                    m = m + 1;
                    Pair(m,1:3) = [SAO(i,1),SAO(j,1),SAO(k,1)];
                    Pair(m,4) = p1;
                    Pair(m,5) = p2;
                    Pair(m,6) = VectorAgl([SAO(i,3),SAO(i,4)],[SAO(j,3),SAO(j,4)],[SAO(k,3),SAO(k,4)]);
                end
            end
        end
    end
end
Pair = Pair(Pair(:,1)>0,:);
save Pair.mat Pair
load SAODATA.mat

ras = SAO(:,3);
decs = SAO(:,4);
mags = SAO(:,3);
%SAO(:,8:10) = [cos(ras).*cos(decs) sin(ras).*cos(decs) sin(decs)];

% Camera Parameter

hozpixel = 1080;
verpixel = 1920;
hoz  = hozpixel/sqrt(hozpixel^2 + verpixel^2);
ver  = verpixel/sqrt(hozpixel^2 + verpixel^2);


hFOV = FOV/2;

SAO(:,8)    =  2*asind(sqrt(sind((SAO(:,4)-ele)/2).^2 + cosd(SAO(:,4)).*cosd(ele).*(sind((SAO(:,3)-azi)/2).^2)));
SAO         = SAO(SAO(:,8)<hFOV,:);

X           = Sphere2Cat(SAO(:,3)',SAO(:,4)');
Y1          = Sphere2Cat(azi,ele);
Y2          = Sphere2Cat(0,0);
% calculate cross and dot products
C = cross(Y1, Y2) ; 
D = dot(Y1, Y2) ;
NP0 = norm(Y1) ; % used for scaling
if ~all(C==0) % check for colinearity    
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    R = sign(D) * (norm(Y2) / NP0) ; % orientation and scaling
end
X1         = (R*X)';

% Plane

X0           = Sphere2Cat(0,0+FOV/2);
N            = Sphere2Cat(0,0);

lend         = size(X1);
for i = 1:lend(1)
    X_new(i,:)       = X1(i,:).*(sum(N.*X0)/sum(N'.*X1(i,:)));
end
Range        = X0.*(sum(N.*X0)/sum(N.*X0))
X_new = X_new./Range(3);
SAO(:,9) = X_new(:,2);
SAO(:,10) = X_new(:,3);


SAO(:,11)   = SAO(:,9)*cosd(ori) - SAO(:,10)*sind(ori);
SAO(:,12)   = SAO(:,9)*sind(ori) + SAO(:,10)*cosd(ori);
plot(SAO(:,11),SAO(:,12),'.g'); hold on
SAO = SAO(abs(SAO(:,11))<ver,:);
SAO = SAO(abs(SAO(:,12))<hoz,:);
plot(SAO(:,11),SAO(:,12),'.b');axis equal;

%grid on
hold on

i = -1:1;
plot(i,hoz*ones(3),'k');hold on
plot(i,-hoz*ones(3),'k');hold on
plot(ver*ones(3),i,'k');hold on
plot(-ver*ones(3),i,'k');hold on
size = size(SAO);
ImgR = zeros(hozpixel,verpixel);
ImgG = zeros(hozpixel,verpixel);
ImgB = zeros(hozpixel,verpixel);
for i = 1:size(1)
     i/size(1)
     Img = generationPointStar(hozpixel,verpixel,((12 - SAO(i,2)))*intensity,SAO(i,12)/hoz,SAO(i,11)/ver);
     ImgR = Img(:,:)*SAO(i,5) + ImgR;
     ImgG = Img(:,:)*SAO(i,6) + ImgG;
     ImgB = Img(:,:)*SAO(i,7) + ImgB;
end
figure
Im(:,:,1) = uint8(ImgR);
Im(:,:,2) = uint8(ImgG);
Im(:,:,3) = uint8(ImgB);

%J = imnoise(Im,'gaussian');
J = Im;
imshow(J)
imwrite(J,'SimSkyNoiseAdded.png');