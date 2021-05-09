%Find true stars
clear all
tic
load Data.mat
%Data = Data(1:6,:);
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
toc
s11 = mode(round(A(:,1)));
s21 = mode(round(A(:,2)*1000))/1000;
s31 = mode(round(A(:,3)*10000))/10000;
s12 = mode(round(A(:,4)));
s22 = mode(round(A(:,5)*1000))/1000;
s32 = mode(round(A(:,6)*10000))/10000;
f(x,y,s1,s2,s3) = s1 + s2*x + s3*y;
Data = Data(abs(Data(:,12)- double(f(Data(:,1),Data(:,2),s11,s21,s31)))<1,:);
toc

sf = fit([Data(:,1),Data(:,2)],Data(:,12),'poly22')
figure
sf2 = fit([Data(:,1),Data(:,2)],Data(:,13),'poly22')
plot(sf,[Data(:,1),Data(:,2)],Data(:,12)); figure
plot(sf2,[Data(:,1),Data(:,2)],Data(:,13))
toc