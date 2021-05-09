function [dis1,dis2,agl] = StarPair(Data,i,j,k)
dis1    = angular_distance(Data(i,8:9),Data(j,8:9));
dis2    = angular_distance(Data(i,8:9),Data(k,8:9));
%angle   = atan2d(norm(cross(vector1,vector2)), dot(vector1,vector2))
%agl     = acosd(dot(vector1,vector2)/(sqrt(vector1(1)^2+vector1(2)^2)*(sqrt(vector2(1)^2+vector2(2)^2))));
agl     = VectorAgl(Data(i,8:9),Data(j,8:9),Data(k,8:9));