function [im,x,y] = generationPointStar(hozpixel,verpixel,std,x,y)
x       = (-x+1)/2*hozpixel;
y       = (-y+1)/2*verpixel;
a       = std*80;
hoz     = 1:hozpixel;
ver     = 1:verpixel;
z1      = (exp(-(hoz - x).^2/(2*std^2)));
z2      = (exp(-(ver - y).^2/(2*std^2)));
im      = (z1'*z2*a);

