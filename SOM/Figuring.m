close all
%% Semi 
figure
plot(KEP(:,1))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,1),'xr')
end
ylabel('Semimajor axis (Meter)')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')
%% Ecc
figure
plot(KEP(:,2))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,2),'xr')
end
ylabel('Eccentricity')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')
%% Inc
figure
plot(KEP(:,3))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,3),'xr')
end
ylabel('Inclination axis (Radian)')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')
%% RAAN
figure
plot(KEP(:,4))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,4),'xr')
end
ylabel('Longitude of the ascending node (Radian)')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')
%%  Argument of periapsis
figure
plot(KEP(:,5))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,5),'xr')
end
ylabel('Argument of periapsis (Radian)')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')

%% Mean anomaly
figure
plot(KEP(:,6))
hold on
for i = 1:10
        plot((Mjd_Epoch(i) - Mjd_Epoch(1))*1440+1,Kep(i,6),'xr')
end
ylabel('Mean anomaly (Radian)')
xlabel('Time since Epoch (minutes)')
legend('Prediction','Reality')