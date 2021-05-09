try 
    load(app.camera_data_dir.Value);
    app.azel_result.Value = '';
    global running
    running = 1;
    filePattern = fullfile('..\Data\StarPos\', '*.mat');
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile('..\Data\StarPos\', baseFileName);
      delete(fullFileName);
    end

    for i = 1:length(Image_Pro.Mjd)
        app.azel_result.Value{i} = sprintf('Calculation for the image #%03d // ', i);     
        Image_pro.UT1   =   Image_Pro.UT1(i);
        Image_pro.AZ    =   Image_Pro.Az(i);
        Image_pro.EL    =   Image_Pro.El(i);
        Image_pro.Exp   =   Image_Pro.ExP(i);
        Image_pro.Nam   =   Image_Pro.Nam(i,:);
        Image_pro.fov   =   Image_Pro.fov(i);
        Image_pro.lon   =   Image_Pro.lon;
        Image_pro.lat   =   Image_Pro.lat;
        [Ra , De , NOS, NOSI]  =   Img2RADE(Image_pro,app,i);
        app.azel_result.Value{end} = [app.azel_result.Value{end} sprintf('[%02d]',NOS)];
        if Ra<(-2*pi)
            AZ(i,2)     =   Ra;
            AZ(i,3)     =   De;
            app.azel_result.Value{end} = [app.azel_result.Value{end} ' Fail in indentifing pattern'];
        else
            AZ(i,2)     =   Ra/180*pi;
            AZ(i,3)     =   De/180*pi;
            app.azel_result.Value{end} = [app.azel_result.Value{end} ' Successfully in indentifing pattern // ' sprintf('AZ: %.5f || EL: %.5f',AZ(i,2),AZ(i,3))];
            try 
                if app.azel_compare.Value == 'On'
                    app.azel_result.Value{end} = [app.azel_result.Value{end} sprintf(' // Different in AZ: %.5f || EL: %.5f',abs(AZ(i,2)- Image_Pro.AZ(i)) ,abs(AZ(i,3)-Image_Pro.EL(i)))];
                end
            catch 
                if app.azel_compare.Value == 'On '
                    app.azel_result.Value{end} = [app.azel_result.Value{end} sprintf(' // Different in AZ: %.5f || EL: %.5f',abs(AZ(i,2)- Image_Pro.AZ(i)) ,abs(AZ(i,3)-Image_Pro.EL(i)))];
                end
            end
        end

        AZ(i,1)         =   (Image_Pro.Mjd(i) + Image_Pro.ExP(i)/86400/2);
        AZ(i,4)         =   Image_Pro.AZ(i);
        AZ(i,5)         =   Image_Pro.EL(i);
        AZ(i,6:11)      =   Image_Pro.Eph(i,:);
        AZ(i,12)        =   Image_Pro.Dis(i);
        AZ(i,13)        =   NOS;
        AZ(i,14)        =   NOSI;

    end
    running = 0;
    fprintf('done\n');
    fprintf('\n');
    save(strcat(app.img_dir.Value,'\AZ.mat'),'AZ');
    app.UseAZELfromAnalysisButton.Enable = 'on';
catch ME
    app.UseAZELfromAnalysisButton.Enable = 'off';
    app.UseAZELfromAnalysisButton.Value = false;
    app.azel_analyze.Value = false;
    app.Lamp.Color = 'red';
    app.Log.Value  = ME.message;
end
