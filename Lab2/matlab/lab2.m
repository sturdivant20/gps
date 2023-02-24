%% Daniel Sturdivant & Andrew Weir | MECH 6970: GPS - LAB 2
clc; clear; close all;

f1 = figure('units','normalized','position',[0.1 0.1 0.3 0.7]);
% f1 = figure('units','normalized','position',[0.1 0.1 0.4 0.5]);
tabs1 = uitabgroup(f1);
tab(1) = uitab(Title='1) SkyPlot');
tab(2) = uitab(Title="2) GeoPlot");
tab(3) = uitab(Title="3) DOP / Num. SV");
tab(4) = uitab(Title="4) ENU Position");
tab(5) = uitab(Title="5) ENU Velocity");
tab(6) = uitab(Title="6) GeoPlot Toomer's");
tab(7) = uitab(Title="7) Dynamic GeoPlot");
tab(8) = uitab(Title="8) Dynamic DOPs");
tab(9) = uitab(Title="9) Speed / Course");
tab(10) = uitab(Title="10) Clock Bias / Drift");
tab(11) = uitab(Title="11) Lab1 SkyPlot");


% constants
load('RCVR_S1_data.mat');
load('RCVR_D1_data.mat');
c = 299792458;  % speed of light

%% PART 1 - STATIC

% toomer' corner
lat0 = 32.606700;
lon0 = -85.481500;
h0 = 215;
% RotMat = [-cosd(lon0)*sind(lat0), -sind(lon0)*sind(lat0), cosd(lat0), 0;
%            sind(lon0)           ,  cosd(lon0)           , 0         , 0;
%            cosd(lon0)*cosd(lat0),  sind(lon0)*cosd(lat0), 0         , 0;
%            0                    ,  0                    , 0         , 1];
RotMat = [-sind(lon0), -sind(lat0)*cosd(lon0), cosd(lat0)*cosd(lon0);
           cosd(lon0), -sind(lat0)*sind(lon0), cosd(lat0)*sind(lon0);
           0         ,  cosd(lat0)           , sind(lon0)];

% unpack sv data and user position
svData1 = svUnpack(RCVR_S1);

az = NaN(30,1819);
el = NaN(30,1819);
time = zeros(1819,1);
lla = zeros(1819,3);
DOP = zeros(4,4,1819);
DOP_ENU = zeros(3,3,1819);
num_sv = zeros(1819,1);
b = zeros(1819,1);
bDot = zeros(1819,1);
x = zeros(1819,3);
xDot = zeros(1819,3);

for i = 1:length(svData1)
    time(i) = svData1{i}.gpsTime - svData1{1}.gpsTime;
    x(i,:) = svData1{i}.x';
    lla(i,:) = svData1{i}.lla;
    xDot(i,:) = svData1{i}.xDot';

    DOP(:,:,i) = svData1{i}.DOP;
    DOP_ENU(:,:,i) = RotMat' * DOP(1:3,1:3,i) * RotMat;

    num_sv(i) = length(svData1{i}.svInUse);
    b(i) = svData1{i}.b;
    bDot(i) = svData1{i}.bDot;

    az(svData1{i}.svInUse,i) = svData1{i}.az;
    el(svData1{i}.svInUse,i) = svData1{i}.el;
end

% SKYPLOT
ax(1) = polaraxes(Parent=tab(1));
title("Skyplot")
labels = ["1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"; "10"; "11"; "12"; "13"; "14"; "15"; "16"; "17"; "18"; "19"; "20"; "21"; "22"; "23"; "24"; "25"; "26"; "27"; "28"; "29"; "30"];
h = skyplot(az(:,1:914)', el(:,1:914)', labels);
% h = skyPlot(az(:,1:914), el(:,1:914), 1:30);
h.LabelFontSize = 16;
h.MaskElevation = 10;

% GEOPLOT
ax(2) = geoaxes(Parent=tab(2));
g = geoplot(ax(2), lla((1:914),1), lla((1:914),2), 'o');
geobasemap(ax(2), "satellite");
geolimits(ax(2), [32.586038, 32.586437], [-85.494653, -85.493987]);

% DOPS PLOT
DOP11 = reshape(DOP(1,1,:), 1, 1819);
DOP22 = reshape(DOP(2,2,:), 1, 1819);
DOP33 = reshape(DOP(3,3,:), 1, 1819);
DOP_ENU11 = reshape(DOP_ENU(1,1,:), 1, 1819);
DOP_ENU22 = reshape(DOP_ENU(2,2,:), 1, 1819);
DOP_ENU33 = reshape(DOP_ENU(3,3,:), 1, 1819);
PDOP = sqrt((DOP11+DOP22+DOP33));
HDOP = sqrt((DOP_ENU11+DOP_ENU22));
VDOP = sqrt((DOP_ENU33));

ax(3) = axes(Parent=tab(3));
sgtitle("Dilution of Precision")
subplot(2,1,2);
hold on;
plot(time(1:914), PDOP(1:914), LineWidth=2.5, DisplayName='PDOP');
plot(time(1:914), HDOP(1:914), LineWidth=2.5, DisplayName='HDOP');
plot(time(1:914), VDOP(1:914), LineWidth=2.5, DisplayName='VDOP');
grid on;
legend;
xlim([0, 914]);
xlabel("Time [s]")
ylabel("DOP")

subplot(2,1,1);
plot(num_sv(1:914), LineWidth=2.5, DisplayName='Num. SV');
grid on;
legend;
xlim([0, 914]);
ylabel("SV in Use")

% ENU PLOT
[xe,xn,xu] = ecef2enu(x(:,1), x(:,2), x(:,3), lat0, lon0, h0, wgs84Ellipsoid("meter"));

ax(4) = axes(Parent=tab(4));
sgtitle("ENU Position Relative to Toomer's Corner")

subplot(3,1,1)
plot(time(1:914), xe(1:914), LineWidth=2.5, DisplayName="East")
grid on;
legend;
xlim([0, 914]);
ylabel("East Position [m]")

subplot(3,1,2)
plot(time(1:914), xn(1:914), LineWidth=2.5, DisplayName="North")
grid on;
legend;
xlim([0, 914]);
ylabel("North Position [m]")

subplot(3,1,3)
plot(time(1:914), xu(1:914), LineWidth=2.5, DisplayName="Up")
grid on;
legend;
xlim([0, 914]);
xlabel("Time [s]")
ylabel("Up Position [m]")

% ENU VELOCITY
[ve,vn,vu] = ecef2enuv(xDot(:,1), xDot(:,2), xDot(:,3), lat0, lon0);

ax(5) = axes(Parent=tab(5));
sgtitle("ENU Velocity")
subplot(3,1,1);
plot(time(1:914), ve(1:914), LineWidth=2.5, DisplayName="East");
grid on;
legend;
xlim([0, 914]);
ylabel("East Velocity [m/s]")

subplot(3,1,2);
plot(time(1:914), vn(1:914), LineWidth=2.5, DisplayName="Norht");
grid on;
legend;
xlim([0, 914]);
ylabel("North Velocity [m/s]")

subplot(3,1,3);
plot(time(1:914), vu(1:914), LineWidth=2.5, DisplayName="Up");
grid on;
legend;
xlim([0, 914]);
ylabel("Up Velocity [m/s]")
xlabel("Time [s]")

% GEOPLOT TOOMER'S
ax(6) = geoaxes(Parent=tab(6));
hold(ax(6), "on");
geoplot(ax(6), lla((1:914),1), lla((1:914),2), 'o', LineWidth=3, MarkerSize=5, DisplayName="Estimates");
geoplot(ax(6), lat0, lon0, '*', LineWidth=3, MarkerSize=10, DisplayName="Toomer's Corner");
geobasemap(ax(6), "streets");


%% PART 2 - DYNAMIC

% unpack sv data and user position
svData2 = svUnpack(RCVR_D1);
idx=~cellfun(@isempty,svData2);
L = sum(idx);

az = NaN(30,L);
el = NaN(30,L);
time = zeros(L,1);
lla = zeros(L,3);
DOP = zeros(4,4,L);
DOP_ENU = zeros(3,3,L);
num_sv = zeros(L,1);
b = zeros(L,1);
bDot = zeros(L,1);
x = zeros(L,3);
xDot = zeros(L,3);

% pull in data
k = 0;
for i = 1:length(idx)

    % only use valid data
    if idx(i) == true
        time(i-k) = svData2{i}.gpsTime - svData2{1}.gpsTime;
        x(i-k,:) = svData2{i}.x';
        lla(i-k,:) = svData2{i}.lla;
        xDot(i-k,:) = svData2{i}.xDot';
    
        DOP(:,:,i-k) = svData2{i}.DOP;
        DOP_ENU(:,:,i-k) = RotMat' * DOP(1:3,1:3,i-k) * RotMat;
    
        num_sv(i-k) = length(svData2{i}.svInUse);
        b(i-k) = svData2{i}.b;
        bDot(i-k) = svData2{i}.bDot;
    
        az(svData2{i}.svInUse,i-k) = svData2{i}.az;
        el(svData2{i}.svInUse,i-k) = svData2{i}.el;
    else
        k = k + 1;
    end

end

% DYNAMIC GEOPLOT
ax(7) = geoaxes(Parent=tab(7));
geoplot(ax(7), lla(:,1), lla(:,2), 'r', LineWidth=2.5);
geobasemap(ax(7), "streets");

% DYNAMIC DOPS
DOP11 = reshape(DOP(1,1,:), 1, L);
DOP22 = reshape(DOP(2,2,:), 1, L);
DOP33 = reshape(DOP(3,3,:), 1, L);
DOP_ENU11 = reshape(DOP_ENU(1,1,:), 1, L);
DOP_ENU22 = reshape(DOP_ENU(2,2,:), 1, L);
DOP_ENU33 = reshape(DOP_ENU(3,3,:), 1, L);
PDOP = sqrt((DOP11+DOP22+DOP33));
HDOP = sqrt((DOP_ENU11+DOP_ENU22));
VDOP = sqrt((DOP_ENU33));

ax(8) = axes(Parent=tab(8));
sgtitle("Dynamic Dilution of Precision")

subplot(2,1,2);
hold on;
plot(time, PDOP, LineWidth=2.5, DisplayName='PDOP');
plot(time, HDOP, LineWidth=2.5, DisplayName='HDOP');
plot(time, VDOP, LineWidth=2.5, DisplayName='VDOP');
grid on;
legend;
ylabel("Error [m]")
xlabel("Time [s]")

subplot(2,1,1);
plot(num_sv, LineWidth=2.5, DisplayName='Num. SV');
grid on;
legend;
ylabel("SV in Use")

% SPEED AND COURSE
[ve,vn,vu] = ecef2enuv(xDot(:,1), xDot(:,2), xDot(:,3), lat0, lon0);
speed = vecnorm([ve,vn,vu],2,2);
course = atan2d(ve, vn);

ax(9) = axes(Parent=tab(9));
sgtitle("Speed and Course")
subplot(2,1,1);
plot(time, speed, LineWidth=2.5, DisplayName="Speed");
grid on;
legend;
ylabel("Velocity Magnitude [m/s]")

subplot(2,1,2);
plot(time, course, LineWidth=2.5, DisplayName="Course");
grid on;
legend;
ylabel("Course [deg]")
xlabel("Time [s]")

% CLOCK BIAS AND DRIFT
ax(10) = axes(Parent=tab(10));
title(ax(10), "Clock Bias and Clock Drift");
subplot(2,1,1);
plot(time, b, LineWidth=2.5, DisplayName="Bias");
grid on;
legend;
ylabel("Bias [m]");

subplot(2,1,2);
plot(time, bDot, LineWidth=2.5, DisplayName="Drift");
grid on;
legend;
xlabel("Time [m]");
ylabel("Drift [m/s]");

% figure
% % plot(time, ve, time, vn, time ,vu)
% plot(time, xDot(:,1), time, xDot(:,2), time ,xDot(:,3))


%% PART 3 - LAB1 EPHEMERIS

ephemerisOut = readMultiGNSSRinex('lab1Data.rnx', 2247);

% unpack sv data 
lat1 = 32.579967;
lon1 = -85.495583;
h1 = 215;
TOW = 2247;
seconds = 24*60*60 + 17*60*60 + 21*60 + 27;

% tmp = cell2mat(struct2cell(ephemerisOut.gps))';
svInUse = [3, 4, 7, 8, 9, 16, 26, 27, 31];

for j = 1:length(svInUse)
    [pos(j,:), vel(j,:), ~] = ...
        svPVTBevly(ephemerisOut.gps(svInUse(j)), seconds, 0.68);
end

[az, el, ~] = ecef2aer(pos(:,1), pos(:,2), pos(:,3), lat1, lon1, h1, wgs84Ellipsoid("meter"));
cond = el > 0;
az(~cond) = NaN;
el(~cond) = NaN;

% LAB1 SKYPLOT
ax(11) = polaraxes(Parent=tab(11));
labels = ["3"; "4"; "7"; "8"; "9"; "16";"26"; "27"; "31"];
h = skyplot(az', el', labels);
% h = skyPlot(az, el, [1:27,29:32]);
h.LabelFontSize = 16;
% h.MaskElevation = 10;

% exportgraphics(tab(1), "../media/static_skyplot.png");
% exportgraphics(tab(2), "../media/static_pos.png");
% exportgraphics(tab(3), "../media/static_dop.png");
% exportgraphics(tab(4), "../media/static_pos_enu.png");
% exportgraphics(tab(5), "../media/static_vel_enu.png");
% exportgraphics(tab(6), "../media/static_geo_toomers.png");
% exportgraphics(tab(7), "../media/dynamic_pos.png");
% exportgraphics(tab(8), "../media/dynamic_dop.png");
% exportgraphics(tab(9), "../media/dynamic_vel_course.png");
% exportgraphics(tab(10), "../media/dynamic_clock.png");
% exportgraphics(tab(11), "../media/lab1_skyplot.png");
