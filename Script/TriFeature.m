%define FOV 
FOV = 40;
hFOV = FOV;

%Load SAO Data
load('SAODATA.mat')
[~, index]= sort(reSAO(:,2),'ascend');
reSAO = reSAO(index,1:4);
length = length(index);


% Pair 3 stars within FOV = 5
m = 0;
Pair = zeros(2895620,6);
for i = 1:length-2
    for j = i+1:length-1
        p1 = angular_distance(reSAO(i,4:-1:3),reSAO(j,4:-1:3));
        if  p1 < hFOV
            for k = j+1:length
                p2 = angular_distance(reSAO(i,4:-1:3),reSAO(k,4:-1:3));
                if p2 <  hFOV
                    m = m + 1;
                    Pair(m,1:3) = [reSAO(i,1),reSAO(j,1),reSAO(k,1)];
                    Pair(m,4) = p1;
                    Pair(m,5) = p2;
                    Pair(m,6) = VectorAgl([reSAO(i,3),reSAO(i,4)],[reSAO(j,3),reSAO(j,4)],[reSAO(k,3),reSAO(k,4)]);
                end
            end
        end
    end
end

Pair = Pair(Pair(:,1)>0,:);
save Pair.mat Pair