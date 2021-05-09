%% Create SAO Database

%% updating Mjd_Epoch
app.Lamp.Enable = 'on';
app.Lamp.Color  = 'blue';
app.Status.Text = 'Generating SAO data ...'; 

try
    addpath ..\Photo_Determination\
    addpath ..\Data\

    fid     = fopen('sao','r');
    YO      = char(app.EpochTime.Value);
    YO      = str2num(YO(7:10))
    i       = 0;
    minStar = app.sao_max.Value;

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
    save(strcat(app.sao_dir.Value,'\SAODATA.mat'),'SAO');
    app.Status.Text = 'Generate SAO data successfully'; 
    app.Lamp.Color = 'green';
    app.Log.Value  = '';
catch ME  
    app.Lamp.Color = 'red';
    app.Status.Text = 'Generate SAO data failed'; 
    app.Log.Value  = ME.message;
end

