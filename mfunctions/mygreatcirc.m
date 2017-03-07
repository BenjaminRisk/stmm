function rng = mygreatcirc(lat1, lon1, lat2, lon2)

% Calculate great circle distance between points on a sphere using the
% Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
% length and has the same units as the radius of the sphere, R.  (If R is
% 1, then RNG is effectively arc length in radians.)

lat1 = lat1*2*pi/360;
lat2 = lat2*2*pi/360;
lon1 = lon1*2*pi/360;
lon2 = lon2*2*pi/360;

a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
 if sum(a > 1)
     warning(['some a are > 1']);
     a(a>1)=1;
 end
if sum(a<0)
    warning(['some a are < 0']);
    a(a<0) = 0;
end

rng = (2 * atan2(sqrt(a),sqrt(1 - a)))*360/2/pi;
