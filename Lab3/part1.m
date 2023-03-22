%% GPS LAB 3: PART 1 - Daniel Sturdivnat & Andrew Weir
clc; close all; clear;
fprintf("<strong>PART 1\n</strong>");

% load("+data/RCVR_S1_data_new2.mat");
% S1 = utils.svUnpack(RCVR_S1, 'sta');
% save("RCVR_S1.mat", "S1");
% 
% load("+data/RCVR_S2_data_new2.mat");
% S2 = utils.svUnpack(RCVR_S2, 'sta');
% save("RCVR_S2.mat", "S2");
% 
% load("+data/RCVR_D1_data_new2.mat");
% D1 = utils.svUnpack(RCVR_D1, 'dyn');
% save("RCVR_D1.mat", "D1");
% 
% load("+data/RCVR_D2_data_new2.mat");
% D2 = utils.svUnpack(RCVR_D2, 'dyn');
% save("RCVR_D2.mat", "D2");

load("RCVR_S1.mat");
c = 299792458;  % speed of light

%% PART A
fprintf("<strong>\n(a)\n</strong>");

% remove unpacked data from struct
% svInUse = [1, 3, 6, 7, 8, 13, 14, 17, 19, 21, 30];
svInUse = [1, 7, 13, 14, 17, 19, 21, 30];
svBool = zeros(32,1);
svBool(svInUse) = 1;

L = length(S1);
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

% psr = psr(svInUse,:);
% dopp = dopp(svInUse,:);
% car = car(svInUse,:);

% plot with satellite view
figure();
geoplot(a_lla(:,1), a_lla(:,2), 'o');


%% PART B
fprintf("<strong>\n(b)\n</strong>");

% carrier smoothed L1 position
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
        [b_x(:,i), b_b(:,i), ~, b_DOP(:,:,i), ~] = ...
            utils.gnssPos(S1{i}.svPos(ismember(svTest, svInUse'),:), p_bar(:,i), 0.5);
        b_lla(i,:,j) = ecef2lla(b_x(:,i)');
    end
end

figure();
geoplot(a_lla(:,1), a_lla(:,2), 'o');
hold on;
geoplot(b_lla(:,1,1), b_lla(:,2,1), 'x');
geoplot(b_lla(:,1,2), b_lla(:,2,2), 'x');
geoplot(b_lla(:,1,3), b_lla(:,2,3), 'x');
legend({'Normal', '2 MIN', '8 MIN', '15 MIN'});
geobasemap satellite;


%% PART C
fprintf("<strong>\n(c)\n</strong>");

% ephemeris ionosphere model
lla0 = [32.5862978360371, -85.4943695259314, 213.567298815586];
lla0_mag = lla0(1) + 0.064*cos(lla0(2) - 1.617);
iono = load("+data/iono_corr_terms_new2.mat");
alpha = [iono.alpha_0; iono.alpha_1; iono.alpha_2; iono.alpha_3];
beta = [iono.beta_0; iono.beta_1; iono.beta_2; iono.beta_3];
AMP = sum(alpha.*lla0_mag);
PER = sum(beta.*lla0_mag);


%% PART D
fprintf("<strong>\n(d)\n</strong>");

% dual frequency
L1 = 1575.42e6;
L2 = 1227.60e6;

psr_IF = (L1^2 / (L1^2 - L2^2))*psr - (L2^2 / (L1^2 - L2^2))*psr2;

d_x = zeros(3,L);
d_b = zeros(L,1);
d_DOP = zeros(4,4,L);
d_lla = zeros(L,3);
for i = 1:L
    svTest = S1{i}.svInUse;
    pp = psr_IF(svTest,i);
    svPos_use = S1{i}.svPos(ismember(svTest, svInUse'),:);
    psr_use = pp(ismember(svTest, svInUse'));
    [d_x(:,i), d_b(:,i), ~, d_DOP(:,:,i), ~] = ...
            utils.gnssPos(svPos_use, psr_use, 0.5);
    d_lla(i,:) = ecef2lla(d_x(:,i)');
end

figure();
geoplot(a_lla(:,1), a_lla(:,2), 'o');
hold on;
geoplot(d_lla(:,1,1), d_lla(:,2,1), 'x');
legend({'Normal', 'IF'});
geobasemap satellite;


%% PART E
fprintf("<strong>\n(e)\n</strong>");

%% PART F
fprintf("<strong>\n(f)\n</strong>");