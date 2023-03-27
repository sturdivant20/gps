%% GPS LAB 3 - PART 4 | Daniel Sturdivant & Andrew Weir
clc; clear; close all;
fprintf("<strong>PART 4\n</strong>");

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Title="A", Parent=tbs);
tab(2) = uitab(Title="A1", Parent=tbs);
tab(3) = uitab(Title="B", Parent=tbs);
tab(4) = uitab(Title="B1", Parent=tbs);
tab(5) = uitab(Title="Comp", Parent=tbs);
tab(6) = uitab(Title="Comp1", Parent=tbs);

load("+data/RCVR_S1.mat");
load("+data/RCVR_D1.mat");

t1 = zeros(1,length(S1));
for i = 1:length(S1)
    t1(i) = S1{i}.gpsTime;
end

t2 = zeros(1,length(D1));
for i = 1:length(D1)
    t2(i) = D1{i}.gpsTime;
end

S1_start = find(t1 == max([t1(1), t2(1)]));
D1_start = find(t2 == max([t1(1), t2(1)]));
S1_finish = find(t1 == min([t1(end), t2(end)]));
D1_finish = find(t2 == min([t1(end), t2(end)]));

for i = 1:(S1_finish-S1_start)
    S1_{i} = S1{i+S1_start-1};
    D1_{i} = D1{i+D1_start-1};
end

S1 = S1_;
D1 = D1_;

clearvars -except S1 D1 tbs tab
L = length(D1);

%% Part A
fprintf("<strong>\n(a)\n</strong>");

% standalone position
x_sta = zeros(3,L);
x_dyn = zeros(3,L);
lla_sta = zeros(3,L);
lla_dyn = zeros(3,L);

for i = 1:L
    % static
    svPos1 = S1{i}.svPos;
    psr1 = S1{i}.L1_psr;
    [x_sta(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos1, psr1, 0.5);
    lla_sta(:,i) = ecef2lla(x_sta(:,i)')';

    % dynamic
    svPos2 = D1{i}.svPos;
    psr2 = D1{i}.L1_psr;
    [x_dyn(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos2, psr2, 0.5);
    lla_dyn(:,i) = ecef2lla(x_dyn(:,i)')';
end

ax = geoaxes(Parent=tab(1));
geoplot(lla_sta(1,:), lla_sta(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn(1,:), lla_dyn(2,:), 'x', LineWidth=2);
legend("Static", "Dynamic");
geobasemap satellite;
title("Standalone Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)

ax = geoaxes(Parent=tab(2));
geoplot(lla_sta(1,:), lla_sta(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn(1,:), lla_dyn(2,:), 'x', LineWidth=2);
legend("Static", "Dynamic");
geobasemap satellite;
geolimits([32.5911, 32.5953], [-85.4864, -85.4803]);
title("Standalone Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% Part B
fprintf("<strong>\n(b)\n</strong>");

% dgps solution
x_sta_dgps = zeros(3,L);
x_dyn_dgps = zeros(3,L);
lla_sta_dgps = zeros(3,L);
lla_dyn_dgps = zeros(3,L);

for i = 1:L
    % static position
    svTest1 = ismember(S1{i}.svInUse, D1{i}.svInUse);
    svPos = S1{i}.svPos;
    psr1 = S1{i}.L1_psr;
    [x_sta_dgps(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos, psr1, 0.5);
    lla_sta_dgps(:,i) = ecef2lla(x_sta_dgps(:,i)')';

    % dynamic dgps position
    svTest2 = ismember(D1{i}.svInUse, S1{i}.svInUse);
    u = svPos(svTest1,:) - x_sta_dgps(:,i)';
    r = sqrt(sum(u.^2,2));

    psr2 = D1{i}.L1_psr(svTest2);
    H = [u./r, ones(size(psr2))];
    y = psr1(svTest1) - psr2;
    dx = pinv(H)*y;
    x_dyn_dgps(:,i) = x_sta_dgps(:,i) + dx(1:3);
    lla_dyn_dgps(:,i) = ecef2lla(x_dyn_dgps(:,i)')';
end

ax = geoaxes(Parent=tab(3));
geoplot(lla_sta_dgps(1,:), lla_sta_dgps(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn_dgps(1,:), lla_dyn_dgps(2,:), 'x', LineWidth=2);
legend("Static", "Dynamic");
geobasemap satellite;
title("DGPS Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)

ax = geoaxes(Parent=tab(4));
geoplot(lla_sta_dgps(1,:), lla_sta_dgps(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn_dgps(1,:), lla_dyn_dgps(2,:), 'x', LineWidth=2);
legend("Static", "Dynamic");
geobasemap satellite;
geolimits([32.5911, 32.5953], [-85.4864, -85.4803]);
title("DGPS Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% plots

% geoplot with both
ax = geoaxes(Parent=tab(5));
geoplot(lla_sta_dgps(1,:), lla_sta_dgps(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn(1,:), lla_dyn(2,:), '^', LineWidth=3);
geoplot(lla_dyn_dgps(1,:), lla_dyn_dgps(2,:), 'x', LineWidth=1.5);
legend("Static", "Dynamic", "Dynamic DGPS", Location="northwest");
geobasemap satellite;
title("Positioning Comparison")
set(findall(gcf,'-property','FontSize'),'FontSize',16)

ax = geoaxes(Parent=tab(6));
geoplot(lla_sta_dgps(1,:), lla_sta_dgps(2,:), 'co', LineWidth=2);
hold on;
geoplot(lla_dyn(1,:), lla_dyn(2,:), '^', LineWidth=3);
geoplot(lla_dyn_dgps(1,:), lla_dyn_dgps(2,:), 'x', LineWidth=1.5);
legend("Static", "Dynamic", "Dynamic DGPS", Location="northwest");
geobasemap satellite;
geolimits([32.5911, 32.5953], [-85.4864, -85.4803]);
title("Positioning Comparison")
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% position difference
f1 = figure(Units='normalized', Position=[1.1, 0.5, 1.3, 0.4]);
tl = tiledlayout(3,1,"TileSpacing","tight", Parent=f1);
xlabel(tl, "Time [s]");
ylabel(tl, "Position Difference [m]");
title(tl, "Position Difference Between Regular and DGPS Position")
ax = axes(Parent=tl);
plot(x_dyn(1,:) - x_dyn_dgps(1,:), 'o', LineWidth=2);
title("X");
grid on;
nexttile;
plot(x_dyn(2,:) - x_dyn_dgps(2,:), 'o', LineWidth=2);
title("Y");
grid on;
nexttile;
plot(x_dyn(3,:) - x_dyn_dgps(3,:), 'o', LineWidth=2);
title("Z");
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%%
exportgraphics(tab(1), "./media/p4_a.png");
exportgraphics(tab(2), "./media/p4_a1.png");
exportgraphics(tab(3), "./media/p4_b.png");
exportgraphics(tab(4), "./media/p4_b1.png");
exportgraphics(f1, "./media/p4_diff.png");
exportgraphics(tab(5), "./media/p4_comp.png");
exportgraphics(tab(6), "./media/p4_comp2.png");
