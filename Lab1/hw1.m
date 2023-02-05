%% Daniel Sturdivant | Andrew Weir | GPS Lab 1
clc; clear; close all;

% Cache ID: GC48DJN
x_cache = [dms2degrees([32,34,0.798*60]), -dms2degrees([85,29,0.735*60])];
x_andrew = [32.579981, -85.495618];
x_daniel = [32.579988, -85.495579];
hgt = distdim(642.58, 'ft', 'm');

% error
x_ned = lla2ned([x_cache,hgt; x_andrew,hgt; x_daniel,hgt], [x_cache,hgt], "flat");
x_rsme = sqrt(sum( (x_ned(1,:)' - x_ned(2:end,:)').^2, 2)) / 2;

fprintf("Cache: %f, %f \nAndrew: %f, %f \nDaniel: %f, %f \n", ...
    x_cache, x_andrew, x_daniel);
fprintf("RSME [m]: %f, %f | %f \n", x_rsme(1:2), norm(x_rsme));

% geoplot
f = figure(Units="normalized",Position=[0 0 1 1]);
ax = geoaxes(FontSize=18);
hold(ax, 'on');
geoplot(ax, x_cache(1), x_cache(2), Color=[0.4660 0.6740 0.1880], ...
    Marker='o', LineWidth=5, MarkerSize=15, LineStyle='none', ...
    DisplayName='Cache');
geoplot(ax, x_andrew(1), x_andrew(2), Color=[0 0.4470 0.7410], ...
    Marker='^', LineWidth=5, MarkerSize=15, LineStyle='none', ...
    DisplayName='Andrew');
geoplot(ax, x_daniel(1), x_daniel(2), Color=[0.6350 0.0780 0.1840], ...
    Marker='square', LineWidth=5, MarkerSize=15, LineStyle='none', ...
    DisplayName='Daniel');
hold(ax, 'off');
geobasemap(ax, 'satellite');
legend(ax);
