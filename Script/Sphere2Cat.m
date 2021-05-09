function y = Sphere2Cat(ra,dec)
y = [cosd(dec).*cosd(ra);cosd(dec).*sind(ra);sind(dec)];
end