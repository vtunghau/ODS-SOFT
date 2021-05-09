function [RA , DE , NOS , NOSI, ME] = Img2RADE(Image_Pro,app,idd)    
    close all
    addpath ..\Photo_Determination\
    load(strcat(app.sao_dir.Value,'\SAODATA.mat'));

    FOVm                = Image_Pro.fov;
    cropsize            = 40;
    
    intensity           = 4;
    
    I                   = imread(Image_Pro.Nam);
    I_gray              = rgb2gray(I);

    I_bin               = im2bw([I_gray;I_gray;I_gray],0.5);
    siz                 = size(I_gray);
    ratio               = siz(2)/siz(1);

    %Find the region
    stats = regionprops('table',I_bin,'Centroid','Area','Eccentricity');
    %filtering data
    stats(stats.Area<intensity,:) = [];
    SAT   = stats(stats.Eccentricity>0.9,:);
    %stats(stats.Area>50,:) = [];
    stats = stats(stats.Eccentricity<0.9,:);
    Centroid = round([stats.Centroid(:,1) stats.Centroid(:,2)]);
    Centroid = [Centroid stats.Area stats.Eccentricity];
    Centroid(Centroid(:,1)<(cropsize/2+1),:) = [];
    Centroid(Centroid(:,2)<(cropsize/2+1),:) = [];
    Centroid(Centroid(:,1)>(siz(2)-cropsize/2-1),:) = [];
    Centroid(Centroid(:,2)>(siz(1)-cropsize/2-1),:) = [];
    Centroid = (double(Centroid));


    
    Img  = zeros(cropsize);
    %Gaussian Fitting
    length_centroid = size(Centroid);
    Data = zeros(length_centroid(1),4);
    for i = 1:length_centroid(1)
        Img = double(I_gray((Centroid(i,2)-cropsize/2):(Centroid(i,2)+cropsize/2-1),(Centroid(i,1)-cropsize/2):(Centroid(i,1)+cropsize/2-1)));
        Data(i,:) = gaussnonlinear2d(Img);
%         subplot(3,4,i)
%         imshow(Img);hold on
%         plot(Data(i,1),Data(i,2),'xr');
%         Str = ['Ecc = ' num2str(Centroid(i,4))];
%         title(Str);
    end

    Data(:,1) = Data(:,1) + Centroid(:,1) - cropsize/2-1;
    Data(:,2) = Data(:,2) + Centroid(:,2) - cropsize/2-1;



    [~,index]= sort(Data(:,4),'descend');
    Data = Data(index,1:4);
%     figure
%     imshow(I); hold on
%     plot(Data(:,1),Data(:,2),'ro')




    Data(:,1)    = Data(:,1)-siz(2)/2;
    Data(:,2)    = Data(:,2)-siz(1)/2;

    FOV          = [1/sqrt(1 + (1/ratio)^2)*FOVm;...
                    1/sqrt(1 + (ratio)^2)*FOVm];

    X            = Sphere2Cat(0,FOVm/2);
    N            = Sphere2Cat(0,0);
    Range        = X.*(sum(N.*X)/sum(N.*X));
    Data(:,5)    = Data(:,1)./(siz(2)/2)*FOV(1)*Range(3)/FOVm;
    Data(:,6)    = Data(:,2)./(siz(1)/2)*FOV(2)*Range(3)/FOVm;
    Data(:,7)    = Range(1);
    [Az, El]     = cart2sph(Data(:,7),Data(:,5),Data(:,6));
    Data(:,8)    = Az/pi*180;
    Data(:,9)    = El/pi*180;
%%  Looping for Pair
    azi          = Image_Pro.AZ/pi*180;
    ele          = Image_Pro.EL/pi*180;
    
    hFOV            = Image_Pro.fov/2;

    minPair         = app.sao_max.Value;
    SAO             = SAO(SAO(:,2)<minPair,:);
    SAO(:,8)        =  2*asind(sqrt(sind((SAO(:,4)-ele)/2).^2 + ...
        cosd(SAO(:,4)).*cosd(ele).*(sind((SAO(:,3)-azi)/2).^2)));
    SAO             = SAO(SAO(:,8)<hFOV,:);
    SAOS            = size(SAO);
    lengt           = SAOS(1);

    STAR            = zeros(7, lengt);

    s = 1;
    Cal      = zeros(29,10e6);
    for i = 1:lengt
        STAR = zeros(7, lengt);
        STAR(6,:) = SAO(i,1);
        STAR(7,:) = SAO(:,1);
        STAR(1,:) = SAO(i,3);
        STAR(2,:) = SAO(i,4);
        STAR(3,:) = SAO(:,3)';
        STAR(4,:) = SAO(:,4)';
        STAR(5,:) =  2*asind(sqrt(sind(abs((STAR(4,:)-STAR(2,:)))/2).^2 + cosd(STAR(4,:)).*cosd(STAR(2,:)).*(sind(abs((STAR(3,:)-STAR(1,:)))/2).^2)));
        STAR = STAR(:,STAR(5,:)<hFOV);
        si   = size(STAR);
        leng = si(2);
        Cal_Array = zeros(20,leng^2);

        for j = 1:leng
            Cal_Array(1,(1 + leng*(j-1)):leng*j) = STAR(6,j);
            Cal_Array(2,(1 + leng*(j-1)):leng*j) = STAR(7,j);
            Cal_Array(3,(1 + leng*(j-1)):leng*j) = STAR(7,:);
            Cal_Array(4,(1 + leng*(j-1)):leng*j) = STAR(5,j);
            Cal_Array(5,(1 + leng*(j-1)):leng*j) = STAR(5,:);
            Cal_Array(6,(1 + leng*(j-1)):leng*j) = STAR(1,j);%az1
            Cal_Array(7,(1 + leng*(j-1)):leng*j) = STAR(3,j);%az2
            Cal_Array(8,(1 + leng*(j-1)):leng*j) = STAR(3,:);
            Cal_Array(9,(1 + leng*(j-1)):leng*j) = STAR(2,j);
            Cal_Array(10,(1 + leng*(j-1)):leng*j) = STAR(4,j);
            Cal_Array(11,(1 + leng*(j-1)):leng*j) = STAR(4,:);
    %         Cal_Array(12:14,(1 + leng*(j-1)):leng*j) ...
    %         = Sphere2Cat(Cal_Array(6,(1 + leng*(j-1)):leng*j),Cal_Array(9,(1 + leng*(j-1)):leng*j));
    %         Cal_Array(15:17,(1 + leng*(j-1)):leng*j) ...
    %         = Sphere2Cat(Cal_Array(7,(1 + leng*(j-1)):leng*j),Cal_Array(10,(1 + leng*(j-1)):leng*j));
    %         Cal_Array(18:20,(1 + leng*(j-1)):leng*j) ...
    %         = Sphere2Cat(Cal_Array(8,(1 + leng*(j-1)):leng*j),Cal_Array(11,(1 + leng*(j-1)):leng*j));
        end
        Cal_Array = Cal_Array(:,Cal_Array(1,:) ~= Cal_Array(2,:));
        Cal_Array = Cal_Array(:,Cal_Array(1,:) ~= Cal_Array(3,:));
        Cal_Array = Cal_Array(:,Cal_Array(2,:) ~= Cal_Array(3,:));
        Cal_Array(12:14,:) = Sphere2Cat(Cal_Array(6,:),Cal_Array(9,:));
        Cal_Array(15:17,:) = Sphere2Cat(Cal_Array(7,:),Cal_Array(10,:));
        Cal_Array(18:20,:) = Sphere2Cat(Cal_Array(8,:),Cal_Array(11,:));
        Cal_Array(21,:) = Cal_Array(15,:) - Cal_Array(12,:);
        Cal_Array(22,:) = Cal_Array(16,:) - Cal_Array(13,:);
        Cal_Array(23,:) = Cal_Array(17,:) - Cal_Array(14,:);
        Cal_Array(24,:) = Cal_Array(18,:) - Cal_Array(12,:);
        Cal_Array(25,:) = Cal_Array(19,:) - Cal_Array(13,:);
        Cal_Array(26,:) = Cal_Array(20,:) - Cal_Array(14,:);
        Cal_Array(4,:) =  2*asind(sqrt(sind(abs(Cal_Array(9,:)-Cal_Array(10,:))/2).^2 + cosd(Cal_Array(9,:)).*cosd(Cal_Array(10,:)).*(sind(abs(Cal_Array(6,:)-Cal_Array(7,:))/2).^2)));
        Cal_Array(5,:) =  2*asind(sqrt(sind(abs(Cal_Array(9,:)-Cal_Array(11,:))/2).^2 + cosd(Cal_Array(9,:)).*cosd(Cal_Array(11,:)).*(sind(abs(Cal_Array(6,:)-Cal_Array(8,:))/2).^2)));
        Cal_Array(28,:) =  2*asind(sqrt(sind(abs(Cal_Array(10,:)-Cal_Array(11,:))/2).^2 + cosd(Cal_Array(10,:)).*cosd(Cal_Array(11,:)).*(sind(abs(Cal_Array(7,:)-Cal_Array(8,:))/2).^2)));
        Nom             = cross(Cal_Array(21:23,:),Cal_Array(24:26,:));
        Cal_Array(27,:) = atan2d(sqrt(Nom(1,:).^2 + Nom(2,:).^2 + Nom(3,:).^2),dot(Cal_Array(21:23,:),Cal_Array(24:26,:)));
        [~,ii] = unique(Cal_Array(27,:));
        Cal_Array = Cal_Array(:,ii);
        %Cal_Array(28,:) =  2*asind(sqrt(sind((Cal_Array(9,:)-Cal_Array(10,:))/2).^2 + cosd(Cal_Array(9,:)).*cosd(Cal_Array(10,:)).*(sind(Cal_Array(6,:)-Cal_Array(7,:))/2).^2));

         Cal_Array(29,:) = Cal_Array(4,:) + Cal_Array(5,:) + Cal_Array(28,:);
        [~,ii] = unique(round(([Cal(29,1:s-1) Cal_Array(29,:)])*10e4)); 
        ii = ii(ii>(s - 1),:);
        Cal_Array = Cal_Array(:,(ii - (s - 1)));
        %Cal_Array = Cal(:,Cal
        e         = size(Cal_Array);
        Cal(:,s:(e(2)+s-1))= Cal_Array;
        s         = s + e(2);
    end


    Cal         = Cal(:,Cal(1,:)~=0);
    Pair(:,1)   = Cal(1,:);
    Pair(:,2)   = Cal(2,:);
    Pair(:,3)   = Cal(3,:);
    Pair(:,4)   = Cal(4,:);
    Pair(:,5)   = Cal(5,:);
    Pair(:,6)   = Cal(27,:);
    Pair(:,7)   = Cal(6,:);
    Pair(:,8)   = Cal(9,:);
    Pair(:,9)   = Cal(7,:);
    Pair(:,10)  = Cal(10,:);
    Pair(:,11)  = Cal(8,:);
    Pair(:,12)  = Cal(11,:);
    
    %% Pivot Search
    load ..\Data\SAODATA.mat
    Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,8)-ele)/2).^2 + cosd(Pair(:,8)).*cosd(ele).*(sind((Pair(:,7)-azi)/2).^2)));
    Pair          = Pair(Pair(:,13)<FOVm/2,:);
    Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,10)-ele)/2).^2 + cosd(Pair(:,10)).*cosd(ele).*(sind((Pair(:,9)-azi)/2).^2)));
    Pair          = Pair(Pair(:,13)<FOVm/2,:);
    Pair(:,13)    =  2*asind(sqrt(sind((Pair(:,12)-ele)/2).^2 + cosd(Pair(:,12)).*cosd(ele).*(sind((Pair(:,11)-azi)/2).^2)));
    Pair          = Pair(Pair(:,13)<FOVm/2,:);
    Pair          = Pair(:,1:6);

    ac_er_ang = 0.1;
    ac_er_dis = 0.01;
    
    DataS = Data;
    fprintf('Number of stars #%d\n',length(Data(:,1)));
    NOS = length(Data(:,1));
    
    try
        for i=1:length(Data(:,1))
            %fprintf('Searching for star #%d\n',i);
            Data(i,10:11) = STAR_NAME(DataS, Pair,ac_er_ang, ac_er_dis,0);
            if Data(i,10) ~= -99999
                DataS = DataS([2:length(Data(:,1)) 1],:);
                Data(i,12:13) = SAO(SAO(:,1)==Data(i,10),3:4);
            end
        end
        Data = Data(Data(:,11)>(max(Data(:,11))/4),:);
    catch ME
        
    end
    
    try
    %% Find True Star
    lendD = size(Data);
    m = 0;
    for i = 1:(lendD-2)
        for j = i+1:(lendD-1)
            for k = j+1:(lendD)
                m = m+1;
                B = [1 Data(i,1) Data(i,2);...
                     1 Data(j,1) Data(j,2);...
                     1 Data(k,1) Data(k,2)];
                C = [Data(i,12);Data(j,12);Data(k,12)];
                C2 = [Data(i,13);Data(j,13);Data(k,13)];
                A(m,1:3) = (B\C);
                A(m,4:6) = (B\C2);
            end
        end
    end
    syms f(x,y,s1,s2,s3)

    s11 = mode(round(A(:,1)));
    s21 = mode(round(A(:,2)*1000))/1000;
    s31 = mode(round(A(:,3)*10000))/10000;
    s12 = mode(round(A(:,4)));
    s22 = mode(round(A(:,5)*1000))/1000;
    s32 = mode(round(A(:,6)*10000))/10000;
    f(x,y,s1,s2,s3) = s1 + s2*x + s3*y;
    Data = Data(abs(Data(:,12)- double(f(Data(:,1),Data(:,2),s11,s21,s31)))<1,:);
    Data = Data(abs(Data(:,13)- double(f(Data(:,1),Data(:,2),s12,s22,s32)))<1,:);
    LenDA= size(Data);
    NOSI = LenDA(1);
    if LenDA(1)>=6
        sf1 = fit([Data(:,1),Data(:,2)],Data(:,12),'poly22');
        sf2 = fit([Data(:,1),Data(:,2)],Data(:,13),'poly22');
        RA  = sf1(0,0);
        DE  = sf2(0,0);
        [RA, DE] = RaDec2AzEl(RA,DE,Image_Pro.lat,Image_Pro.lon, Image_Pro.UT1);

%         figure
%         plot(sf1,[Data(:,1),Data(:,2)],Data(:,12)); 
%         figure
%         plot(sf2,[Data(:,1),Data(:,2)],Data(:,13))
    else
        fprintf('Rejected case\n');
        RA  = -999;
        DE  = -999;
    end
    
    Data(:,1) = Data(:,1)+app.width.Value/2;
    Data(:,2) = Data(:,2)+app.height.Value/2;
%     plot(Data(:,1),Data(:,2),'ob')
%     str = (Data(:,10));
%     text(Data(:,1)+20,Data(:,2),num2str(str),'Color','w');
    save(['..\Data\StarPos\' sprintf('StarPos%04d.mat', idd)],'Data'); 
    catch ME
        NOSI = 0;
        fprintf('Rejected case\n');
        RA  = -999;
        DE  = -999;
    end
end