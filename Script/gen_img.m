function [outdata, Ps] = gen_img(Azel_Data,DT1,DT2,app)
    addpath     ..\Photo_Determination\
    load(strcat(app.sao_dir.Value,'\SAODATA.mat'));
    FOV                 = app.fov.Value;
    intensity           = app.intensity.Value ;
    intensity_max       = app.sao_max.Value + 3;
    hozpixel            = app.height.Value;
    verpixel            = app.width.Value;
    azi                 = (Azel_Data(DT2,11)/pi*180 +  Azel_Data(DT1,11)/pi*180)/2;
    ele                 = (Azel_Data(DT2,12)/pi*180 +  Azel_Data(DT1,12)/pi*180)/2;
    ori                 = app.ori.Value;
    % Camera Parameter

    
    hoz  = hozpixel/sqrt(hozpixel^2 + verpixel^2);
    ver  = verpixel/sqrt(hozpixel^2 + verpixel^2);


    hFOV = FOV/2;

    SAO(:,8)    =  2*asind(sqrt(sind((SAO(:,4)-ele)/2).^2 + cosd(SAO(:,4)).*cosd(ele).*(sind((SAO(:,3)-azi)/2).^2)));
    SAO         = SAO(SAO(:,8)<hFOV,:);

    X           = Sphere2Cat(SAO(:,3)',SAO(:,4)');
    Y1          = Sphere2Cat(azi,ele);
    Y2          = Sphere2Cat(0,0);

    % calculate cross and dot products
    C           = cross(Y1, Y2) ; 
    D           = dot(Y1, Y2) ;
    NP0         = norm(Y1) ; % used for scaling
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
    Range        = X0.*(sum(N.*X0)/sum(N.*X0));
    X_new = X_new./Range(3);
    SAO(:,9) = X_new(:,2);
    SAO(:,10) = X_new(:,3);


    SAO(:,11)   = (SAO(:,9)*cosd(ori) - SAO(:,10)*sind(ori));
    SAO(:,12)   = SAO(:,9)*sind(ori) + SAO(:,10)*cosd(ori);
%     plot(SAO(:,11),SAO(:,12),'.g'); hold on
    SAO = SAO(abs(SAO(:,11))<ver,:);
    SAO = SAO(abs(SAO(:,12))<hoz,:);
%     plot(-SAO(:,11),SAO(:,12),'.b','MarkerSize',15);
    
%     close all
%     figure
%     for i = 1:length(SAO(:,11))
%         plot(-SAO(i,11),SAO(i,12),'.','MarkerEdgeColor',[SAO(i,5) SAO(i,6) SAO(i,7)],'MarkerSize',((intensity_max - SAO(i,2)))*4);
%         hold on
%     axis equal
%     end
%     set(gca,'Color', [150 150 150]/255)
%     xlim([-ver ver]);
%     ylim([-hoz hoz]);
    
    
    
    SAT                 = [Azel_Data(DT1,11)/pi*180 Azel_Data(DT1,12)/pi*180;
                           Azel_Data(DT2,11)/pi*180 Azel_Data(DT2,12)/pi*180];
    XSAT                = Sphere2Cat(SAT(:,1)',SAT(:,2)');
    XSAT_1              = (R*XSAT)';
    XSAT_new(1,:)       = XSAT_1(1,:).*(sum(N.*X0)/sum(N'.*XSAT_1(1,:)));
    XSAT_new(2,:)       = XSAT_1(2,:).*(sum(N.*X0)/sum(N'.*XSAT_1(2,:)));
    XSAT_new            = XSAT_new./Range(3);
    XSAT_new(:,2)       = XSAT_new(:,2)*cosd(ori) - XSAT_new(:,2)*sind(ori);
    XSAT_new(:,3)       = XSAT_new(:,3)*sind(ori) + XSAT_new(:,3)*cosd(ori);
    out                 = path2im(XSAT_new(1,2:3), XSAT_new(2,2:3), verpixel, hozpixel)*200;

%     
%     hold off
%     plot(-XSAT_new(:,2),XSAT_new(:,3),'Color',[1 1 1],'LineWidth',2)
%     hold on
%     set(gca,'Color', [150 150 150]/255)
%     axis equal
%     xlim([-ver ver]);
%     ylim([-hoz hoz]);
    
    
    
%     grid on
%     hold on

    i = -1:1;
%     plot(i,hoz*ones(3),'k');hold on
%     plot(i,-hoz*ones(3),'k');hold on
%     plot(ver*ones(3),i,'k');hold on
%     plot(-ver*ones(3),i,'k');hold on
    sizeSAO = size(SAO);

    ImgR = zeros(hozpixel,verpixel);
    ImgG = zeros(hozpixel,verpixel);
    ImgB = zeros(hozpixel,verpixel);
    for i = 1:sizeSAO(1)
         [Img,x,y] = generationPointStar(hozpixel,verpixel,((intensity_max - SAO(i,2)))*intensity,SAO(i,12)/hoz,SAO(i,11)/ver);
         Ps(i,:)   = [x,y];
         ImgR = Img(:,:)*SAO(i,5) + ImgR;
         ImgG = Img(:,:)*SAO(i,6) + ImgG;
         ImgB = Img(:,:)*SAO(i,7) + ImgB;
    end

%     figure
%     imshow(uint8(out))
    
    noise = imnoise(uint8(out),'gaussian');

%     noise = uint8(zeros(hozpixel,verpixel));
    Im(:,:,1) = uint8(ImgR) + noise;
    Im(:,:,2) = uint8(ImgG) + noise;
    Im(:,:,3) = uint8(ImgB) + noise;
%     figure
%     imshow(Im)
   
    outdata = Im;
end