%% GPS LAB 3: PART 1 - Daniel Sturdivnat & Andrew Weir
clc; close all; clear;
fprintf("<strong>PART 2\n</strong>");

f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab("A", Parent=tbs);
tab(2) = uitab("B", Parent=tbs);
tab(3) = uitab("C", Parent=tbs);
tab(4) = uitab("D", Parent=tbs);
tab(5) = uitab("F", Parent=tbs);

c = 299792458;  % speed of light
L1 = 1575.42e6;
L2 = 1227.60e6;
lambda1 = c / L1;
lambda2 = c / L2;

svInUse = [1, 7, 13, 14, 17, 19, 21, 30]';
% svInUse = [1, 7, 14, 17, 30]';
T = length(svInUse);
svInUse_ = zeros(32,1);
svInUse_(svInUse) = 1;

load("+data/RCVR_S1.mat");
r1.L = length(S1);

r1.gpsTime = zeros(r1.L,1);
r1.L1psr = nan(32,r1.L);
r1.L1psr_var = nan(32,r1.L);
r1.L1dopp = nan(32,r1.L);
r1.L1car = nan(32,r1.L);
r1.L1car_var = nan(32,r1.L);
r1.L2psr = nan(32,r1.L);
r1.L2psr_var = nan(32,r1.L);
r1.L2dopp = nan(32,r1.L);
r1.L2car = nan(32,r1.L);
r1.L2car_var = nan(32,r1.L);

for i = 1:r1.L
    r1.gpsTime(i) = S1{i}.gpsTime;
    svTest = S1{i}.svInUse;

    r1.L1psr(svInUse,i) = S1{i}.L1_psr(ismember(svTest, svInUse));
    r1.L1psr_var(svInUse,i) = S1{i}.L1_psr_var(ismember(svTest, svInUse));
    r1.L1dopp(svInUse,i) = S1{i}.L1_dopp(ismember(svTest, svInUse));
    r1.L1car(svInUse,i) = S1{i}.L1_car(ismember(svTest, svInUse));
    r1.L1car_var(svInUse,i) = S1{i}.L1_car_var(ismember(svTest, svInUse));

    r1.L2psr(svInUse,i) = S1{i}.L2_psr(ismember(svTest, svInUse));
    r1.L2psr_var(svInUse,i) = S1{i}.L2_psr_var(ismember(svTest, svInUse));
    r1.L2dopp(svInUse,i) = S1{i}.L2_dopp(ismember(svTest, svInUse));
    r1.L2car(svInUse,i) = S1{i}.L2_car(ismember(svTest, svInUse));
    r1.L2car_var(svInUse,i) = S1{i}.L2_car_var(ismember(svTest, svInUse));
end

load("+data/RCVR_S2.mat");
r2.L = length(S2);

r2.gpsTime = zeros(r2.L,1);
r2.L1psr = nan(32,r2.L);
r2.L1psr_var = nan(32,r2.L);
r2.L1dopp = nan(32,r2.L);
r2.L1car = nan(32,r2.L);
r2.L1car_var = nan(32,r1.L);
r2.L2psr = nan(32,r2.L);
r2.L2psr_var = nan(32,r2.L);
r2.L2dopp = nan(32,r2.L);
r2.L2car = nan(32,r1.L);
r2.L2car_var = nan(32,r1.L);

for i = 1:r2.L
    r2.gpsTime(i) = S2{i}.gpsTime;
    svTest = S2{i}.svInUse;

    r2.L1psr(svInUse,i) = S2{i}.L1_psr(ismember(svTest, svInUse));
    r2.L1psr_var(svInUse,i) = S2{i}.L1_psr_var(ismember(svTest, svInUse));
    r2.L1dopp(svInUse,i) = S2{i}.L1_dopp(ismember(svTest, svInUse));
    r2.L1car(svInUse,i) = S2{i}.L1_car(ismember(svTest, svInUse));
    r2.L1car_var(svInUse,i) = S2{i}.L1_car_var(ismember(svTest, svInUse));

    r2.L2psr(svInUse,i) = S2{i}.L2_psr(ismember(svTest, svInUse));
    r2.L2psr_var(svInUse,i) = S2{i}.L2_psr_var(ismember(svTest, svInUse));
    r2.L2dopp(svInUse,i) = S2{i}.L2_dopp(ismember(svTest, svInUse));
    r2.L2car(svInUse,i) = S2{i}.L2_car(ismember(svTest, svInUse));
    r2.L2car_var(svInUse,i) = S2{i}.L2_car_var(ismember(svTest, svInUse));
end


% only use portions where times are the same
r1.Start = find(r1.gpsTime == max([r1.gpsTime(1), r2.gpsTime(1)]));
r2.Start = find(r2.gpsTime == max([r1.gpsTime(1), r2.gpsTime(1)]));
r1.Finish = find(r1.gpsTime == min([r1.gpsTime(end), r2.gpsTime(end)]));
r2.Finish = find(r2.gpsTime == min([r1.gpsTime(end), r2.gpsTime(end)]));

r1.L = r1.Finish - r1.Start + 1;
r1.gpsTime = r1.gpsTime(r1.Start:r1.Finish);
r1.L1psr = r1.L1psr(:, r1.Start:r1.Finish);
r1.L1psr_var = r1.L1psr_var(:, r1.Start:r1.Finish);
r1.L1dopp = r1.L1dopp(:, r1.Start:r1.Finish);
r1.L1car = r1.L1car(:, r1.Start:r1.Finish);
r1.L2psr = r1.L2psr(:, r1.Start:r1.Finish);
r1.L2psr_var = r1.L2psr_var(:, r1.Start:r1.Finish);
r1.L2dopp = r1.L2dopp(:, r1.Start:r1.Finish);
r1.L2car = r1.L2car(:, r1.Start:r1.Finish);
r1.L2car_var = r1.L2car_var(:, r1.Start:r1.Finish);

r2.L = r2.Finish - r2.Start + 1;
r2.gpsTime = r2.gpsTime(r2.Start:r2.Finish);
r2.L1psr = r2.L1psr(:, r2.Start:r2.Finish);
r2.L1psr_var = r2.L1psr_var(:, r2.Start:r2.Finish);
r2.L1dopp = r2.L1dopp(:, r2.Start:r2.Finish);
r2.L1car = r2.L1car(:, r2.Start:r2.Finish);
r2.L2psr = r2.L2psr(:, r2.Start:r2.Finish);
r2.L2psr_var = r2.L2psr_var(:, r2.Start:r2.Finish);
r2.L2dopp = r2.L2dopp(:, r2.Start:r2.Finish);
r2.L2car = r2.L2car(:, r2.Start:r2.Finish);
r2.L2car_var = r2.L2car_var(:, r2.Start:r2.Finish);


%% PART A
fprintf("<strong>\n(a)\n</strong>");

% path of the two receivers
a_x1 = zeros(3,r1.L);
a_lla1 = zeros(3,r1.L);
a_x2 = zeros(3,r2.L);
a_lla2 = zeros(3,r2.L);
a_dist = zeros(1, r1.L);
for i = 1:r1.L % same ar r2.L
    svTest1 = S1{i+r1.Start-1}.svInUse;
    svPos1 = S1{i+r1.Start-1}.svPos(ismember(svTest1, svInUse),:);

    svTest2 = S2{i+r2.Start-1}.svInUse;
    svPos2 = S2{i+r2.Start-1}.svPos(ismember(svTest2, svInUse),:);

    [a_x1(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos1, r1.L1psr(svInUse,i), 0.5);
    a_lla1(:,i) = ecef2lla(a_x1(:,i)')';

    [a_x2(:,i), ~, ~, ~, ~] = utils.gnssPos(svPos2, r2.L1psr(svInUse,i), 0.5);
    a_lla2(:,i) = ecef2lla(a_x2(:,i)')';

    a_dist(i) = norm(a_x1(:,i) - a_x2(:,i));
end

a_mu = mean(a_dist);
a_std = std(a_dist);
fprintf("std = %f\n", a_std);
fprintf("mean = %f\n", a_mu);

ax = geoaxes(Parent=tab(1));
geoplot(a_lla1(1,:), a_lla1(2,:), 'o');
hold on;
geoplot(a_lla2(1,:), a_lla2(2,:), 'x');
legend("RCVR1", "RCVR2");
geobasemap satellite;
title("Regular Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART B
fprintf("<strong>\n(b)\n</strong>");

% psuedorange DGPS (ref = r2, user = r1)
b_x = zeros(3,r1.L);
b_lla = zeros(3,r1.L);
b_dist = zeros(1,r1.L);

for i = 1:r1.L
    % unit vector and range from ref (known)
    svTest2 = S2{i+r2.Start-1}.svInUse;
    svPos2 = S2{i+r2.Start-1}.svPos(ismember(svTest2, svInUse),:);
    u = svPos2 - a_x2(:,i)';
    r = sqrt(sum(u.^2,2));

    % dgps
    dPsuedo = r1.L1psr(svInUse,i) - r2.L1psr(svInUse,i);
    H = [u./r, ones(T,1)];
    dx = pinv(H)*dPsuedo;
    b_x(:,i) = a_x2(:,i) + dx(1:3);
    b_lla(:,i) = ecef2lla(b_x(:,i)')';

    % distance
    b_dist(i) = norm(a_x2(:,i) - b_x(:,i));
end

b_mu = mean(b_dist);
b_std = std(b_dist);
fprintf("std = %f\n", b_std);
fprintf("mean = %f\n", b_mu);

ax = geoaxes(Parent=tab(2));
geoplot(a_lla2(1,:), a_lla2(2,:), 'o');
hold on;
geoplot(b_lla(1,:), b_lla(2,:), 'x');
legend("RCVR2", "RCVR1 DGPS");
geobasemap satellite;
title("DGPS Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART C
fprintf("<strong>\n(c)\n</strong>");

% carrier smoothed DGPS
M = [2,8,15].*60;
c_x = zeros(3,r1.L,3);
c_lla = zeros(3,r1.L,3);
c_dist = zeros(3,r1.L);

for j = 1:3
    dCar = zeros(T,r1.L);
    dPsuedo = zeros(T,r1.L);

    for i = 1:r1.L
        % unit vector and range from ref (known)
        svTest2 = S2{i+r2.Start-1}.svInUse;
        svPos2 = S2{i+r2.Start-1}.svPos(ismember(svTest2, svInUse),:);
        u = svPos2 - a_x2(:,i)';
        r = sqrt(sum(u.^2,2));
    
        % smoothing
        dP = r2.L1psr(svInUse,i) - r1.L1psr(svInUse,i);
        dCar(:,i) = lambda1.*abs(r2.L1car(svInUse,i)) - lambda1.*abs(r1.L1car(svInUse,i));
        if i == 1
            dPsuedo = r2.L1psr(svInUse,i) - r1.L1psr(svInUse,i);
        elseif i <= M(j)
            dPsuedo(:,i) = (1/i)*dP + ((i-1)/i)*(dPsuedo(:,i-1) + dCar(:,i) - dCar(:,i-1));
        else
            dPsuedo(:,i) = (1/M(j))*dP + ((M(j)-1)/M(j))*(dPsuedo(:,i-1) + dCar(:,i) - dCar(:,i-1));
        end

        % dgps
        H = [u./r, ones(T,1)];
        dx = pinv(H)*dPsuedo(:,i);
        c_x(:,i,j) = a_x2(:,i) + dx(1:3);
        c_lla(:,i,j) = ecef2lla(c_x(:,i,j)')';
    
        % distance
        c_dist(j,i) = norm(a_x2(:,i) - c_x(:,i,j));
    end
end

c_mu = mean(c_dist,2);
c_std = std(c_dist,[],2);
fprintf("std = %f, %f, %f\n", c_std);
fprintf("mean = %f, %f, %f\n", c_mu);

ax = geoaxes(Parent=tab(3));
geoplot(a_lla2(1,:), a_lla2(2,:), 'o', LineWidth=5);
hold on;
geoplot(c_lla(1,:,1), c_lla(2,:,1), 'x', LineWidth=3);
geoplot(c_lla(1,:,2), c_lla(2,:,2), 'x', LineWidth=1.5);
geoplot(c_lla(1,:,3), c_lla(2,:,3), 'x', LineWidth=0.5);
legend("RCVR2", "RCVR1 2 MIN" , "RCVR1 8 MIN", "RCVR1 15 MIN");
geobasemap satellite;
title("Carrier Smoothed DGPS Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART D
fprintf("<strong>\n(d)\n</strong>");

addpath(genpath("Lambda"));

% carrier (rtk) DGPS
d_x = zeros(3,r1.L);
d_lla = zeros(3,r1.L);
d_dist = zeros(1,r1.L);
d_N = zeros(2*T,r1.L);
d_N2 = zeros(2*T,r1.L);
d_N_cov = zeros(2*T,2*T,r1.L);

for i = 1:r1.L
    svTest2 = S2{i+r2.Start-1}.svInUse;
    svPos2 = S2{i+r2.Start-1}.svPos(ismember(svTest2, svInUse),:);
    u = svPos2 - a_x2(:,i)';
    r = sqrt(sum(u.^2,2));
    
    % recursive least squares for N
    R = diag([r1.L1psr_var(svInUse,i) + r2.L1psr_var(svInUse,i); ...
              r1.L1car_var(svInUse,i) + r2.L1car_var(svInUse,i); ...
              r1.L2car_var(svInUse,i) + r2.L2car_var(svInUse,i)]);

    m = [r2.L1psr(svInUse,i) - r1.L1psr(svInUse,i); ...
         abs(r2.L1car(svInUse,i).*lambda1) - abs(r1.L1car(svInUse,i).*lambda1); ...
         abs(r2.L2car(svInUse,i).*lambda2) - abs(r1.L2car(svInUse,i).*lambda2)];
    
    G = [u./r, ones(T,1); ...
         u./r, ones(T,1); ...
         u./r, ones(T,1)];

    LAM = [zeros(T,T)              , zeros(T,T); ...
           lambda1.*diag(ones(T,1)), zeros(T,T); ...
           zeros(T,T)              , lambda2.*diag(ones(T,1))];

    L = null(G')';

%     d_N(:,i) = pinv(L*LAM)*(L*m);
%     d_cov_N(:,:,i) = ((L*LAM)'*(L*Q^-1*L')*(L*LAM))^-1;

    RR = L*R^-1*L';
    HH = L*LAM;
    yy = L*m;
    d_N(:,i) = (HH'*RR*HH)^-1*HH'*RR*yy;
    d_N_cov(:,:,i) = (HH'*RR*HH)^-1;
%     if i == 1
%         Q = HH' * RR * HH;
%         d_N(:,1) = (Q^-1 * HH' * RR * yy);
% %         d_N_cov(:,:,1) = Q^-1;
% %         d_N(:,1) = zeros(2*T,1);
%         d_N_cov(:,:,1) = eye(length(d_N(:,1))).*10;
%     else
%         Q = Q + (HH' * RR * HH);
%         dx = Q^-1 * HH' * (yy - HH*d_N(:,i-1));
% 
%         d_N(:,i) = d_N(:,i-1) + dx;
%         d_N_cov(:,:,i) = Q^-1;
%     end

    % lambda method
    [d_N2,~,~,~,~,~,~] = LAMBDA(d_N(:,i), d_N_cov(:,:,i));

    % find position
    y = abs(r2.L1car(svInUse,i)*lambda1) - abs(r1.L1car(svInUse,i)*lambda1) - (lambda1.*round(d_N2(1:T,1)));
    H = [u./r, ones(T,1)];
    dx = pinv(H)*y;
    d_x(:,i) = a_x2(:,i) + dx(1:3);
    d_lla(:,i) = ecef2lla(d_x(:,i)')';

    d_dist(i) = norm(a_x2(:,i) - d_x(:,i));
end

ax = geoaxes(Parent=tab(4));
geoplot(a_lla2(1,:), a_lla2(2,:), 'o');
hold on;
geoplot(d_lla(1,:), d_lla(2,:), 'x');
legend("RCVR2", "RCVR1 RTK DGPS");
geobasemap satellite;
title("RTK DGPS Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)

d_mu = mean(d_dist,2);
d_std = std(d_dist,[],2);
fprintf("std = %f\n", d_std);
fprintf("mean = %f\n", d_mu);


%% distance plot

ax = axes(Parent=tab(5));
hold on;
plot(1:1987, a_dist, 'o', LineWidth=2.5);
plot(1:1987, b_dist, 'o', LineWidth=1.5);
plot(1:1987, c_dist(1,:), 'x', LineWidth=2);
plot(1:1987, c_dist(2,:), 'x', LineWidth=2);
plot(1:1987, c_dist(3,:), 'x', LineWidth=2);
plot(1:1987, d_dist(1,:), 'square', LineWidth=2);
legend("Regular", "Code DGPS", "Carrier Smoothed DGPS 2 MIN", "Carrier Smoothed DGPS 8 MIN", "Carrier Smoothed DGPS 15 MIN", "RTK DGPS")
title("Static Receiver Positioning Difference")
xlabel("Time [s]");
ylabel("Position Difference [m]");
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% 

% exportgraphics(tab(1), "./media/p2_a.png");
% exportgraphics(tab(2), "./media/p2_b.png");
% exportgraphics(tab(3), "./media/p2_c.png");
% exportgraphics(tab(4), "./media/p2_d.png");
% exportgraphics(tab(5), "./media/p2_comp.png");
