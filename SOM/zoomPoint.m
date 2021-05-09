function out = zoomPoint(Y0,size)
xlim([Y0(1)-size Y0(1)+size]);
ylim([Y0(2)-size Y0(2)+size]);
zlim([Y0(3)-size Y0(3)+size]);
end