function geom_lat = geog2geom(geog_lat, geog_long, height, decimalYear)
%INPUTS:
%geog_lat: Geographic latitude, specified as a scalar or vector, in degrees. North latitude is positive and south latitude is negative.
%geog_long: Geographic longitude specified as a scalar or vector, in degrees. East longitude is positive and west longitude is negative.
%height: Distance from the surface of the Earth, specified as a scalar or vector, in meters.
%decimalYear: Year, in decimal format, specified as a scalar or vector. This value can have any fraction of the year that has already passed.
%The geog_lat, geog_long, height, and decimalYear arguments must all be the same size (scalars or vectors).
%OUTPUT
%geom_lat: Geomagnetic latitude, specified as a scalar or vector, in degrees. North latitude is positive and south latitude is negative.
%The size of geom_lat is same as the size of each of the inputs
%number of data points
n=length(geog_lat);
geom_lat=nan(size(geog_lat));
for i=1:n
    %First compute the magnetic dip I(also known as the magnetic inclination or dip angle)
    [~,~,~,I,~,~,~,~,~,~] = igrfmagm(height(i),geog_lat(i),geog_long(i),decimalYear(i));
    %Then compute the magnetic latitude (geom_lat) as: tand(MLAT)=0.5*tand(I)
    geom_lat(i)=atand(0.5*tand(I));
end