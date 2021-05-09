Test_loop
% Searching Name of PiPivotStarot Star in Frame
addpath ..\Photo_Determination\
load ..\Data\Data.mat
load ..\Data\dataPair7.mat
load ..\Data\SAODATA.mat
load ..\Data\Azel_Data.mat
load ..\Data\Input.mat
FOVm                = Input.fov;
%load Pair.mat
%load SAODATA.mat

azi = Azel_Data(3,2)/pi*180
ele = Azel_Data(3,3)/pi*180
Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,8)-ele)/2).^2 + cosd(Pair(:,8)).*cosd(ele).*(sind((Pair(:,7)-azi)/2).^2)));
Pair         = Pair(Pair(:,13)<FOVm/2,:);
Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,10)-ele)/2).^2 + cosd(Pair(:,10)).*cosd(ele).*(sind((Pair(:,9)-azi)/2).^2)));
Pair         = Pair(Pair(:,13)<FOVm/2,:);
Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,12)-ele)/2).^2 + cosd(Pair(:,12)).*cosd(ele).*(sind((Pair(:,11)-azi)/2).^2)));
Pair         = Pair(Pair(:,13)<FOVm/2,:);
Pair = Pair(:,1:6);

ac_er_ang = 0.02;
ac_er_dis = 0.01;



tic
%% Identificate The Pivot Star in The pattern Image
%Data(1,3:4) = STAR_NAME(Data, Pair,ac_er_ang, ac_er_dis,0);
%Data(i,5:6) = reSAO(reSAO(:,1)==Data(i,3),3:4);
DataS = Data;

for i=1:length(Data(:,1))
    i
    Data(i,10:11) = STAR_NAME(DataS, Pair,ac_er_ang, ac_er_dis,0);
    DataS = DataS([2:length(Data(:,1)) 1],:);
    Data(i,12:13) = SAO(SAO(:,1)==Data(i,10),3:4);
end
Data = Data(Data(:,11)>(max(Data(:,11))/3),:);
save Data.mat Data

toc






