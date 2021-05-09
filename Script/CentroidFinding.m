clc;
close all;
clear all;


addpath ..\Photo_Determination\
addpath ..\Data\
load    ..\Data\Input.mat\
load    ..\Data\Ps.mat
%2020
tic
FOVm                = Input.fov;
cropsize            = 16;
ratio               = 16/9;
intensity           = 1;

I = imread('..\Image\Img.png');

%I = imread('SimSkyNoiseAdded.png');
I_gray = rgb2gray(I);
I_bin  = im2bw(I,0.6);
size = size(I_gray);


%Find the region
stats = regionprops('table',I_bin,'Centroid','Area','Eccentricity');
%filtering data
stats(stats.Area<intensity,:) = [];
SAT   = stats(stats.Eccentricity>0.9,:);
%stats(stats.Area>50,:) = [];
stats = stats(stats.Eccentricity<0.9,:);
Centroid = round([stats.Centroid(:,1) stats.Centroid(:,2)]);
Centroid(Centroid(:,1)<cropsize/2,:) = [];
Centroid(Centroid(:,2)<cropsize/2,:) = [];
Centroid(Centroid(:,1)>(size(2)-cropsize/2),:) = [];
Centroid(Centroid(:,2)>(size(1)-cropsize/2),:) = [];
Centroid = double(Centroid);


Data = zeros(length(Centroid),4);
Img = zeros(cropsize);
%Gaussian Fitting
for i = 1:length(Centroid)  
    Img = double(I_gray((Centroid(i,2)-cropsize/2):(Centroid(i,2)+cropsize/2-1),...
                        (Centroid(i,1)-cropsize/2):(Centroid(i,1)+cropsize/2-1)));
    Data(i,:) = gaussnonlinear2d(Img);  
end

Data(:,1) = Data(:,1) + Centroid(:,1) - cropsize/2-1;
Data(:,2) = Data(:,2) + Centroid(:,2) - cropsize/2-1;



[~,index]= sort(Data(:,4),'descend');
Data = Data(index,1:4);
imshow(I); hold on
%save Data.mat Data
plot(Data(:,1),Data(:,2),'or')



Ps2(:,1)     = Data(:,1);
Ps2(:,2)     = Data(:,2);
Data(:,1) = Data(:,1)-size(2)/2;
Data(:,2) = Data(:,2)-size(1)/2;

FOV  = [1/sqrt(1 + (1/ratio)^2)*FOVm;...
        1/sqrt(1 + (ratio)^2)*FOVm];

X            = Sphere2Cat(0,FOVm/2);
N            = Sphere2Cat(0,0);
Range        = X.*(sum(N.*X)/sum(N.*X));
Data(:,5)    = Data(:,1)./(size(2)/2)Range(3)/FOVm;
Data(:,6)    = Data(:,2)./(size(1)/2)Range(3)/FOVm;

Data(:,7)    = Range(1);
[Az, El]     = cart2sph(Data(:,7),Data(:,5),Data(:,6));
Data(:,8)    = Az/pi*180;
Data(:,9)    = El/pi*180;
save ..\Data\Data.mat Data



run('PivotSearch.m')
run('FindTrueStar.m')

Data(:,1) = Data(:,1)+1920/2;
Data(:,2) = Data(:,2)+1080/2;
str = (Data(:,10));
text(Data(:,1)+20,Data(:,2),num2str(str),'Color','w');

toc