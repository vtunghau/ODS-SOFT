%% updating Mjd_Epoch
year = app.year.Value;
mon  = app.month.Value;
day  = app.day.Value;
hr   = app.hour.Value;
minute = app.min.Value;
sec  = app.sec.Value;


Mjd_Epoch = Mjday(year,mon,day,hr,minute,sec);
%Mjd_Epoch = fix(Mjd_Epoch) + mod(epoch,1);
app.Mjd_epoch.Value = Mjd_Epoch;
app.Status.Text = 'Updating time successfully'; 
