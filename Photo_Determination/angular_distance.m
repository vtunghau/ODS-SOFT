function angle = angular_distance(point1,point2)
    point1 = point1/180*pi;
    point2 = point2/180*pi;
    angle = 2*asin(sqrt(sin(abs((point1(1)-point2(1)))/2)^2 + cos(point1(1))*...
           cos(point2(1))*sin(abs(point1(2)-point2(2))/2)^2))/pi*180;