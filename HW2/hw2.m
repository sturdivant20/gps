%% Daniel Sturdivant | GPS HW2
clc; clear; close all;

f = figure('units','normalized','position',[0.1 0.1 0.3 0.7]);
tabs = uitabgroup(f);
tab(1) = uitab(Title='1) y=a');
tab(2) = uitab(Title="2) LS Fitting");
tab(3) = uitab(Title="3) Ranging");
tab(4) = uitab(Title="4) Psuedoranging");
tab(5) = uitab(Title="5) Elevation Mask");
tab(6) = uitab(Title="5) Skyplot");


%% Problem 1
fprintf("PROBLEM 1 \n")

N = 1000;   % num monte carlo runs
M = 500;    % max num meas
sig_exp = sqrt(1 ./ (1:M));  % expected variance
sig = zeros(M,1);
for m = 1:M
    err = zeros(m*N,1);
    for n = 1:m:m*N
        err(n:n+m-1) = ones(m,1)-randn(m,1);
    end
    sig(m) = sqrt(var(err) / m); % calculate variance
end

ax = axes(Parent=tab(1));
hold(ax, "on");
plot(ax, 1:M, sig_exp, Color="#0072BD", MarkerSize=8, LineWidth=3, DisplayName='Expected');
plot(ax, 1:M, sig, ':', Color="#EDB120", MarkerSize=8, LineWidth=3, DisplayName='Monte Carlo');
grid(ax, "on");
title(ax, "\textbf{Monte Carlo Standard Deviation Estimates of y=1 $\sim\mathcal{N}$(0,1)}", Interpreter="latex");
xlabel(ax, "Number of Measurements");
ylabel(ax, "Standard Deviation (\sigma)");
legend(ax);
ax.FontSize = 16;
ax.FontName = 'Arial';

exportgraphics(tab(1), './media/p1.png');
fprintf("A bunch of ploots\n");


%% Problem 2
fprintf("\nPROBLEM 2 \n");

x = [0    ; 1    ; 2    ; 3    ; 4    ];
y = [0.181; 2.680; 3.467; 3.101; 3.437];
sigma = 0.4;

H{1} = [ones([5,1]), x];                % y = a + bx
H{2} = [ones([5,1]), x, x.^2];          % y = a + bx + cx^2
H{3} = [ones([5,1]), x, x.^2, x.^3];    % y = a + bx + cx^2 + dx^3
X{1} = inv(H{1}'*H{1})*H{1}'*y;
X{2} = inv(H{2}'*H{2})*H{2}'*y;
X{3} = inv(H{3}'*H{3})*H{3}'*y;
P{1} = sigma^2 * inv(H{1}'*H{1});
P{2} = sigma^2 * inv(H{2}'*H{2});
P{3} = sigma^2 * inv(H{3}'*H{3});

xx = linspace(0,4,1000);
Y{1} = X{1}(1) + X{1}(2)*xx;
Y{2} = X{2}(1) + X{2}(2)*xx + X{2}(3)*xx.^2;
Y{3} = X{3}(1) + X{3}(2)*xx + X{3}(3)*xx.^2 + X{3}(4)*xx.^3;

ax = axes(Parent=tab(2));
hold(ax, "on");
plot(x,y, 'o-', Color="#77AC30", LineWidth=4, MarkerSize=12, DisplayName="Reference Data");
plot(xx,Y{1}, '--', Color="#0072BD",LineWidth=3, DisplayName="1st order");
plot(xx,Y{2}, ':', Color="#A2142F", LineWidth=3, DisplayName="2nd order");
plot(xx,Y{3}, '-.', Color="#EDB120", LineWidth=3, DisplayName="3rd order");
grid(ax, "on");
legend(ax, Location="southeast");
title("\textbf{Least Squares Fitting}", Interpreter="latex")
xlabel("x");
ylabel("y");
ax.FontSize = 16;

exportgraphics(tab(2), './media/p2.png');
fprintf("'a' for each estimate: %f | %f | %f \n", X{1}(1), X{2}(1), X{3}(1));


%% Problem 3
fprintf("\nPROBLEM 3 \n");

a  = [0 ; 10; 0 ; 10];
b  = [0 ; 0 ; 10; 10];
r2 = [25; 65; 45; 85];
sigma = 0.5;

% perfect measurements
error = Inf;
x = [1;1];
i = 0;
while error > 1e-3
    i = i+1;
    H = -2 .* [a-x(1), b-x(2)];             % linear measurement model
    y = r2 - ((a-x(1)).^2 + (b-x(2)).^2);   % measurement estimate
%     H = [-(a-x(1)), -(b-x(2))] ./ sqrt((a-x(1)).^2 + (b-x(2)).^2);
%     y = sqrt(r2) - sqrt((a-x(1)).^2 + (b-x(2)).^2);
    dx = inv(H'*H)*H'*y;                    % LS
    x = x + dx;
    error = norm(dx);
end
P = sigma^2 .* inv(H'*H);

fprintf("Num Iter: %d \n", i);
fprintf("[X,Y] = [%f, %f] \n", x);
% fprintf("P = [%f, %f; \n    %f, %f] \n", P);

% monte carlo version
x_hat = zeros(2,1000);
P_hat = zeros(2,2,1000);
j = 0;
for i = 1:10000
    error = inf;
    x = [1;1];
    noise = sigma .* randn(size(r2)); % range measurement noise
    while error > 1e-3
        j = j+1;
        H = -2 .* [a-x(1), b-x(2)];
        y = (r2+noise) - ((a-x(1)).^2 + (b-x(2)).^2);
%         H = [-(a-x(1)), -(b-x(2))] ./ sqrt((a-x(1)).^2 + (b-x(2)).^2);
%         y = (sqrt(r2)+noise) - sqrt((a-x(1)).^2 + (b-x(2)).^2);
        dx = inv(H'*H)*H'*y;
        x = x + dx;
        error = norm(dx);
    end
    x_hat(:,i) = x;
    P_hat(:,:,i) = sigma^2 .* inv(H'*H);
end
sig_P = mean(P_hat,3);
% sig_MC = cov((x_hat)');
sig_MC = std(x_hat');

fprintf("sig_exp(x), sig_exp(y) = %f, %f \n", sqrt(P(1)), sqrt(P(4)));
% fprintf("sig(x), sig(y) = %f, %f \n", sqrt(sig_MC(1)), sqrt(sig_MC(4)));
fprintf("sig(x), sig(y) = %f, %f \n", sig_MC);


%% Problem 4
fprintf("\nPROBLEM 4 \n");

% read text file
opts = detectImportOptions('sv_pos_one_epoch.txt');
opts = setvartype(opts,"double");
opts.Delimiter = ',';

opts.DataLines = [3, 11];
gps_data = readmatrix('sv_pos_one_epoch.txt',opts);
opts.DataLines = [14, 19];
soop_data = readmatrix('sv_pos_one_epoch.txt', opts);

sigma = 0.5;
ax = geoaxes(Parent=tab(4));
hold(ax, "on");
geobasemap(ax, 'satellite');

x = gnssPsuedo(gps_data(1:9, 2:4), gps_data(1:9, 5), sigma);
lla = ecef2lla(x');
RotMat = [-cosd(lla(2))*sind(lla(1)), -sind(lla(2))*sind(lla(1)), cosd(lla(1)), 0;
           sind(lla(2))             ,  cosd(lla(2))             , 0           , 0;
           cosd(lla(2))*cosd(lla(1)),  sind(lla(2))*cosd(lla(1)), 0           , 0;
           0                        ,  0                        , 0           , 1];

% sv 1-4
x_sv = gps_data(1:4, 2:4);
rho = gps_data(1:4, 5);
[x,b,P,DOP,i] = gnssPsuedo(x_sv, rho, sigma);
fprintf("Using 4 satellites: (%d iterations) \n", i);
fprintf("x = [%f, %f, %f] \n  = [%f, %f, %f] \nb = %f \n", x, ecef2lla(x'), b);
lla = ecef2lla(x');
geoplot(lla(1), lla(2), 'o', Color="#77AC30", MarkerSize=5, LineWidth=2, DisplayName='4 SV');
Q = RotMat * DOP * RotMat';
PDOP = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));
fprintf("PDOP = %f \nHDOP = %f \nVDOP = %f \n", PDOP, HDOP, VDOP);

% sv 1-4 no clock bias
[x,b,P,DOP,i] = gnssPsuedo(x_sv, rho-100, sigma);
fprintf("Using 4 satellites with corrected/no clock bias: (%d iterations) \n", i);
fprintf("x = [%f, %f, %f] \n  = [%f, %f, %f] \nb = %f \n", x, ecef2lla(x'), b);
lla = ecef2lla(x');
geoplot(lla(1), lla(2), 'o', Color="#0072BD", MarkerSize=5, LineWidth=2, DisplayName='4 SV, b=0');
Q = RotMat * DOP * RotMat';
PDOP = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));
fprintf("PDOP = %f \nHDOP = %f \nVDOP = %f \n", PDOP, HDOP, VDOP);

% sv 1-9
x_sv = gps_data(:, 2:4);
rho = gps_data(:, 5);
[x,b,P,DOP,i] = gnssPsuedo(x_sv, rho, sigma);
fprintf("Using all 9 satellites: (%d iterations) \n", i);
fprintf("x = [%f, %f, %f] \n  = [%f, %f, %f] \nb = %f \n", x, ecef2lla(x'), b);
lla = ecef2lla(x');
geoplot(lla(1), lla(2), 'o', Color="#EDB120", MarkerSize=5, LineWidth=2, DisplayName='9 SV');
Q = RotMat * DOP * RotMat';
PDOP = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));
fprintf("PDOP = %f \nHDOP = %f \nVDOP = %f \n", PDOP, HDOP, VDOP);

% 2 gnss and 2 soop
x_sv = [gps_data(6:7, 2:4); soop_data(2:3, 2:4)];
rho = [gps_data(6:7, 5); soop_data(2:3, 5)];
[x,b,P,DOP,i] = gnssPsuedo(x_sv, rho, sigma);
% fprintf("Using 2 GNSS and 2 SOOP satellites: (%d iterations) \n", i);
fprintf("x = [%f, %f, %f] \n  = [%f, %f, %f] \nb = %f \n", x, ecef2lla(x'), b);
lla_soop = ecef2lla(x');
geoplot(lla_soop(1), lla_soop(2), 'o', MarkerSize=5, LineWidth=2, DisplayName='2 SOOP');
Q = RotMat * DOP * RotMat';
PDOP = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));
fprintf("PDOP = %f \nHDOP = %f \nVDOP = %f \n", PDOP, HDOP, VDOP);

% sv 1-9 wiht initial condition
x_sv = gps_data(:, 2:4);
rho = gps_data(:, 5);
[x,b,P,DOP,i] = gnssPsuedo(x_sv, rho, sigma, [423000; -5362000; 3417000; 0]);
fprintf("Using all 9 satellites with IC: (%d iterations) \n", i);
fprintf("x = [%f, %f, %f] \n  = [%f, %f, %f] \nb = %f \n", x, ecef2lla(x'), b);
lla = ecef2lla(x');
hold on;
geoplot(lla(1), lla(2), 'o', Color="#A2142F", MarkerSize=5, LineWidth=2, DisplayName='9 SV w/ IC');
Q = RotMat * DOP * RotMat';
PDOP = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));
fprintf("PDOP = %f \nHDOP = %f \nVDOP = %f \n", PDOP, HDOP, VDOP);

legend(ax);

%% Problem 5
fprintf("\nPROBLEM 5 \n");

A1 = 5*10^-9;

% known reciever position
sigma = 0.5;
x_sv = gps_data(:, 2:4);
rho = gps_data(:, 5);
[x_known,b_known,P_known,~,~] = gnssPsuedo(x_sv, rho, sigma);

RotMat = [-cosd(lla(2))*sind(lla(1)), -sind(lla(2))*sind(lla(1)), cosd(lla(1)), 0;
           sind(lla(2))             ,  cosd(lla(2))             , 0           , 0;
           cosd(lla(2))*cosd(lla(1)),  sind(lla(2))*cosd(lla(1)), 0           , 0;
           0                        ,  0                        , 0           , 1];

% find azimuth and elevation of each satellite
d_xyz = x_sv' - x_known;
% az = atan2d( d_xyz(2,:) , d_xyz(1,:) )';
% el = atand( -d_xyz(3,:) ./ sqrt(d_xyz(1,:).^2 + d_xyz(2,:).^2) )';
lla_known = ecef2lla(x_known');
[az, el, R] = ecef2aer(x_sv(:,1), x_sv(:,2), x_sv(:,3), lla_known(1), lla_known(2), lla_known(3), wgs84Ellipsoid("meter"));

% loop from horizontal to vertical
for theta = 0:90
    
    i = theta + 1;

    % check elevation angles
    sv_use = x_sv(abs(el) > theta, :);
    rho_use = rho(abs(el) > theta);
    num_sv(i) = sum(abs(el) > theta);

    if num_sv(i) < 4
        max_el = theta-1;
        num_sv = num_sv(1:end-1);
%         [x(:,i),b(i),P(:,:,i),itr(i)] = gnssPsuedo(sv_use, rho_use, sigma);
%         I(i) = A1 *H (1 + 16*(0.53 - theta/180)^3);
        break;
    end

    % ionoshpere error
    I(i) = A1 * (1 + 16*(0.53 - (theta/180))^3);

    % calculate location
    sigma = el(abs(el) > theta).^(-1);
    [x(:,i),b(i),P(:,:,i),DOP(:,:,i),itr(i)] = gnssPsuedo(sv_use, rho_use, sigma);
    Q(:,:,i) = RotMat * DOP(:,:,i) * RotMat';

end

theta = (1:length(num_sv)) - 1;
axes(Parent=tab(5));
subplot(3,1,1);
plot(theta, num_sv, LineWidth=3, Color="#77AC30");
grid on;
title("\textbf{Number of Satellites in Use}", Interpreter="latex");
ylabel("Number of Satellites");
ax = gca;
ax.FontSize = 16;
subplot(3,1,2);
hold on;
P11 = reshape(DOP(1,1,:), 1, length(num_sv));
P22 = reshape(DOP(2,2,:), 1, length(num_sv));
P33 = reshape(DOP(3,3,:), 1, length(num_sv));
Q11 = reshape(Q(1,1,:), 1, length(num_sv));
Q22 = reshape(Q(2,2,:), 1, length(num_sv));
Q33 = reshape(Q(3,3,:), 1, length(num_sv));
PDOP = sqrt((P11+P22+P33));
HDOP = sqrt((Q11+Q22));
VDOP = sqrt((Q33));
% plot(theta, P11, LineWidth=2, Color="#0072BD", DisplayName='X Error');
% plot(theta, P22, LineWidth=2, Color="#EDB120", DisplayName='Y Error');
% plot(theta, P33, LineWidth=2, Color="#A2142F", DisplayName='Z Error');
plot(theta, PDOP, LineWidth=2, Color="#4DBEEE", DisplayName='PDOP');
plot(theta, HDOP, LineWidth=2, Color="#EDB120", DisplayName='HDOP');
plot(theta, VDOP, LineWidth=2, Color="#A2142F", DisplayName='VDOP');
title("\textbf{Dilution of Precision}", Interpreter="latex");
ylabel("DOPs [m]");
ylim([0,45]);
yticks([0,15,30,45]);
grid on;
ax = gca;
ax.FontSize = 16;
legend(Location='northwest');
subplot(3,1,3);
plot(theta, I, LineWidth=3, Color="#7E2F8E");
title("\textbf{Ionosphere Error}", Interpreter="latex");
xlabel("Elevation [deg]");
ylabel("Delay [s]");
grid on;
ax = gca;
ax.FontSize = 16;

exportgraphics(tab(5), './media/p5.png');

ax = polaraxes(Parent=tab(6));
skyplot(az, el, 'MaskElevation',max_el);