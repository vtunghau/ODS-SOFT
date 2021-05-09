app.Lamp.Color  = 'blue';
app.Status.Text = 'Generating Image . . .';

try 

    load(strcat(app.img_dir.Value,'\Azel_Data.mat'));
    load(strcat(app.img_dir.Value,'\Sta.mat'));
    app.img_item.Items     = {};
    app.img_item.Enable    = 'off';
    for i = 1:app.noo.Value
        k = i*3;
        fprintf('Generated image #%d\n',i);
        outdata            = gen_img(Azel_Data,k-2,k,app);
        Num                = sprintf('%04d',i);
        name_file          = [app.img_dir.Value '\Image\' Num '.png'];
        imwrite(outdata, name_file);
        Image_Pro.Mjd(i)   = Azel_Data(k-1,1);
        Image_Pro.AZ(i)    = Azel_Data(k-1,2);
        Image_Pro.EL(i)    = Azel_Data(k-1,3);
        Image_Pro.Az(i)    = Azel_Data(k-1,11);
        Image_Pro.El(i)    = Azel_Data(k-1,12);
        Image_Pro.Dis(i)   = Azel_Data(k-1,10);
        Image_Pro.ExP(i)   = app.obs_time.Value;
        Image_Pro.Nam(i,:) = name_file;
        Image_Pro.fov(i)   = app.fov.Value;
        Image_Pro.Eph(i,:) = Azel_Data(k-1,4:9);
        Image_Pro.UT1(i)   = Azel_Data(k-1,13);
        Image_Pro.lon      = Sta.lon;
        Image_Pro.lat      = Sta.lat;
        app.img_item.Items{end + 1}   = strcat('Image #',Num);
        app.Log.Value{end} = sprintf('Generated image #%d',i);
        app.Log.Value{end+1} = '';
%         scroll(app.Log,length(app.Log.Value));
    end
    save(strcat(app.img_dir.Value,'\Image_Pro.mat'),'Image_Pro');
    app.img_item.Enable     = 'on';
    app.Lamp.Color          = 'green';
    app.Status.Text         = 'Generated Image successfully';

catch ME
    app.Log.Value   = ME.message;
    app.Lamp.Color  = 'red';
    app.Status.Text = 'Generated Image failed';
    
end
