%% GPS LAB 3: PART 1 - Daniel Sturdivnat & Andrew Weir
clc; close all; clear;
fprintf("<strong>PART 1\n</strong>");

f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab("A", Parent=tbs);
tab(2) = uitab("B", Parent=tbs);
tab(3) = uitab("C", Parent=tbs);
tab(4) = uitab("D", Parent=tbs);
tab(5) = uitab("F", Parent=tbs);

% load("+data/RCVR_S1_data_new2.mat");
% S1 = utils.svUnpack(RCVR_S1, 'sta');
% save("+data/RCVR_S1.mat", "S1");
% 
% load("+data/RCVR_S2_data_new2.mat");
% S2 = utils.svUnpack(RCVR_S2, 'sta');
% save("+data/RCVR_S2.mat", "S2");
% 
% load("+data/RCVR_D1_data_new2.mat");
% D1 = utils.svUnpack(RCVR_D1, 'dyn');
% save("+data/RCVR_D1.mat", "D1");
% 
% load("+data/RCVR_D2_data_new2.mat");
% D2 = utils.svUnpack(RCVR_D2, 'dyn');
% save("+data/RCVR_D2.mat", "D2");

load("+data/RCVR_S1.mat");

c = 299792458;  % speed of light
L1 = 1575.42e6;
L2 = 1227.60e6;
lambda1 = c / L1;
lambda2 = c / L2;

L = length(S1);
% svInUse = [1, 3, 6, 7, 8, 13, 14, 17, 19, 21, 30];
svInUse = [1, 7, 13, 14, 17, 19, 21, 30];

%% PART A
fprintf("<strong>\n(a)\n</strong>");

% remove unpacked data from struct
svBool = zeros(32,1);
svBool(svInUse) = 1;

gpsTime = zeros(L,1);
num_sv = zeros(L,1);
b = zeros(L,1);
x = zeros(3,L);
a_lla = zeros(L,3);
DOP = zeros(4,4,L);
psr = nan(32,L);
dopp = nan(32,L);
car = nan(32,L);
psr2 = nan(32,L);
dopp2 = nan(32,L);

for i = 1:L
    gpsTime(i) = S1{i}.gpsTime;
    num_sv(i) = length(S1{i}.svInUse);
    psr(S1{i}.svInUse,i) = S1{i}.L1_psr;
    dopp(S1{i}.svInUse,i) = S1{i}.L1_dopp;
    car(S1{i}.svInUse,i) = S1{i}.L1_car;
    psr2(S1{i}.svInUse,i) = S1{i}.L2_psr;
    dopp2(S1{i}.svInUse,i) = S1{i}.L2_dopp;

    svTest = S1{i}.svInUse; 
    % generate L1 position
    [a_x(:,i), a_b(:,i), ~, a_DOP(:,:,i), ~] = ...
        utils.gnssPos(S1{i}.svPos(ismember(svTest, svInUse),:), psr(svInUse,i), 0.5);
    a_lla(i,:) = ecef2lla(a_x(:,i)');
end

a_std = std(a_x, [], 2);
a_mu = mean(a_lla, 1);
fprintf("Mean [LLA] = %f, %f, %f\n", a_mu);
fprintf("Std [XYZ] = %f, %f, %f\n", a_std);

% psr = psr(svInUse,:);
% dopp = dopp(svInUse,:);
% car = car(svInUse,:);

% plot with satellite view
ax = geoaxes(Parent=tab(1));
geoplot(a_lla(:,1), a_lla(:,2), 'o');
geobasemap satellite;
title("Regular Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART B
fprintf("<strong>\n(b)\n</strong>");

% carrier smoothed L1 position
b_x = zeros(3,L,3);
b_b = zeros(L,1);
b_DOP = zeros(4,4,L);
b_lla = zeros(L,3,3);

p_bar = nan(8,1999);
M = [2,8,15].*60;
for j = 1:3
    for i = 1:L
        if i == 1
            p_bar(:,1) = psr(svInUse,1);
        elseif i < M(j)
            p_bar(:,i) = (1/i)*psr(svInUse,i) + ...
                ((i-1)/i)*(p_bar(:,i-1) + (car(svInUse,i) - car(svInUse,i-1))*-(c/1575.42e6));
        else
            p_bar(:,i) = (1/M(j))*psr(svInUse,i) + ...
                ((M(j)-1)/M(j))*(p_bar(:,i-1) + (car(svInUse,i) - car(svInUse,i-1))*-(c/1575.42e6));
        end
    
        svTest = S1{i}.svInUse;
        % generate L1 position
        [b_x(:,i,j), b_b(:,i), ~, b_DOP(:,:,i), ~] = ...
            utils.gnssPos(S1{i}.svPos(ismember(svTest, svInUse'),:), p_bar(:,i), 0.5);
        b_lla(i,:,j) = ecef2lla(b_x(:,i,j)');
    end
end

b_std1 = std(b_x(:,:,1), [], 2);
b_mu1 = mean(b_lla(:,:,1), 1);
b_std2 = std(b_x(:,:,2), [], 2);
b_mu2 = mean(b_lla(:,:,2), 1);
b_std3 = std(b_x(:,:,3), [], 2);
b_mu3 = mean(b_lla(:,:,3), 1);
fprintf("Mean [LLA] = %f, %f, %f\n", b_mu1);
fprintf("Std [XYZ] = %f, %f, %f\n", b_std1);
fprintf("Mean [LLA] = %f, %f, %f\n", b_mu2);
fprintf("Std [XYZ] = %f, %f, %f\n", b_std2);
fprintf("Mean [LLA] = %f, %f, %f\n", b_mu3);
fprintf("Std [XYZ] = %f, %f, %f\n", b_std3);

ax = geoaxes(Parent=tab(2));
geoplot(a_lla(:,1), a_lla(:,2), 'o');
hold on;
geoplot(b_lla(:,1,1), b_lla(:,2,1), 'x');
geoplot(b_lla(:,1,2), b_lla(:,2,2), 'x');
geoplot(b_lla(:,1,3), b_lla(:,2,3), 'x');
legend({'Normal', '2 MIN', '8 MIN', '15 MIN'});
geobasemap satellite;
title("Carrier Smoothed Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART C
fprintf("<strong>\n(c)\n</strong>");

% ephemeris ionosphere model
lla0 = [32.5862978360371, -85.4943695259314, 213.567298815586];
iono = load("+data/iono_corr_terms_new2.mat");
alpha = [iono.alpha_0; iono.alpha_1; iono.alpha_2; iono.alpha_3];
beta = [iono.beta_0; iono.beta_1; iono.beta_2; iono.beta_3];
c_x = zeros(3,L);
c_b = zeros(L,1);
c_DOP = zeros(4,4,L);
c_lla = zeros(L,3);
T_iono = zeros(8,L);

for i = 1:L
    % 1) PULL DATA FROM SV EPHEMERIS
    sv = ismember(S1{i}.svInUse', svInUse);
    gpsTime = S1{i}.gpsTime;
    psuedo = S1{i}.L1_psr(sv);
    svPos = S1{i}.svPos(sv,:);

    [ecef0, ~, ~, ~, ~] = utils.gnssPos(svPos, psuedo, 0.5);
    lla0 = ecef2lla(ecef0');

    [A,E,~] = ecef2aer(svPos(:,1), svPos(:,2), svPos(:,3), lla0(1), lla0(2), lla0(3), wgs84Ellipsoid("meter"));
    A = deg2rad(A);
    E = deg2rad(E);

    % 2) IONOSPERE MODEL
    psi = (0.0137 / (E + 0.11)) - 0.022;
    
    phi_i = deg2rad(lla0(1)) + psi*cos(A);
    if phi_i > 0.416
        phi = 0.416;
    elseif phi_i < -0.416
        phi_i = -0.416;
    end
    
    lambda_i = deg2rad(lla0(2)) + (psi*sin(A) / cos(phi_i));
    phi_m = phi_i + 0.064*cos(lambda_i - 1.617);
    
    t = 4.32e4 * lambda_i + gpsTime;
    if t > 86400
        t = t - 86400;
    elseif t < 0
        t = t + 86400;
    end
    
    F = 1 + 16*(0.53 - E).^3;
    
    PER = sum(beta.*phi_m);
    if PER < 72000
        PER = 72000;
    end
    
    AMP = sum(alpha.*phi_m);
    if AMP < 0
        AMP = 0;
    end
    
    x = 2*pi*(t - 50400) / PER;
    if abs(x) < 1.57
        T_iono(:,i) = F * (5e-9 + AMP*(1 - (x^2/2) + (x^4/24)));
    else
        T_iono(:,i) = F * 5e-9;
    end

    % 3) CORRECT PSUEDORANGE AND RECALCULATE
    psuedo = psuedo - c.*T_iono(:,i);
    [c_x(:,i), c_b(:,i), ~, c_DOP(:,:,i), ~] = utils.gnssPos(svPos, psuedo, 0.5);
    c_lla(i,:) = ecef2lla(c_x(:,i)');
end

c_std = std(c_x, [], 2);
c_mu = mean(c_lla, 1);
fprintf("Mean [LLA] = %f, %f, %f\n", c_mu);
fprintf("Std [XYZ] = %f, %f, %f\n", c_std);

ax = geoaxes(Parent=tab(3));
geoplot(a_lla(:,1), a_lla(:,2), 'o');
hold on;
geoplot(c_lla(:,1), c_lla(:,2), 'x');
legend({'Normal', 'Iono. Model'});
geobasemap satellite;
title("Ionospheric Model Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART D
fprintf("<strong>\n(d)\n</strong>");

% dual frequency
L1 = 1575.42e6;
L2 = 1227.60e6;

psr_IF = (L1^2 / (L1^2 - L2^2)).*psr - (L2^2 / (L1^2 - L2^2)).*psr2;

d_x = zeros(3,L);
d_b = zeros(L,1);
d_DOP = zeros(4,4,L);
d_lla = zeros(L,3);
for i = 1:L
    psuedo = psr_IF(S1{i}.svInUse,i);
    svPos = S1{i}.svPos(ismember(S1{i}.svInUse, svInUse'),:);
    psuedo = psuedo(ismember(S1{i}.svInUse, svInUse'));
    [d_x(:,i), d_b(:,i), ~, d_DOP(:,:,i), ~] = utils.gnssPos(svPos, psuedo, 0.5);
    d_lla(i,:) = ecef2lla(d_x(:,i)');
end

d_std = std(d_x, [], 2);
d_mu = mean(d_lla, 1);
fprintf("Mean [LLA] = %f, %f, %f\n", d_mu);
fprintf("Std [XYZ] = %f, %f, %f\n", d_std);

ax = geoaxes(Parent=tab(4));
geoplot(a_lla(:,1), a_lla(:,2), 'o');
hold on;
geoplot(d_lla(:,1), d_lla(:,2), 'x');
legend({'Normal', 'Dual Freq.'});
geobasemap satellite;
title("Dual Frequency Positioning")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% PART E
fprintf("<strong>\n(e)\n</strong>");

% NOAA better ephemeris

%% Part F
fprintf("<strong>\n(f)\n</strong>");

% all tests
ax = geoaxes(Parent=tab(5));
geoplot(a_lla(:,1), a_lla(:,2), 'o', LineWidth=9);
hold on;
geoplot(b_lla(:,1,1), b_lla(:,2,1), '^', LineWidth=1.5);
geoplot(b_lla(:,1,2), b_lla(:,2,2), '^', LineWidth=1.5);
geoplot(b_lla(:,1,3), b_lla(:,2,3), '^', LineWidth=1.5);
geoplot(c_lla(:,1), c_lla(:,2), 'square', LineWidth=1.5);
geoplot(d_lla(:,1), d_lla(:,2), 'x', LineWidth=2.5);
geobasemap satellite;
legend("Regular", "Carrier Smoothed 2 MIN", "Carrier Smoothed 8 MIN", "Carrier Smoothed 15 MIN", "Iono. Model", "Dual Freq.")
title("Positioning Comparison")
set(findall(gcf,'-property','FontSize'),'FontSize',16)


%% 

% exportgraphics(tab(1), "./media/p1_a.png");
% exportgraphics(tab(2), "./media/p1_b.png");
% exportgraphics(tab(3), "./media/p1_c.png");
% exportgraphics(tab(4), "./media/p1_d.png");
% exportgraphics(tab(5), "./media/p1_f.png");
