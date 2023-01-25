%% Daniel Sturdivant | MECH 6970: GPS | HW1
clc; close all; clear;

%% PART 1
fprintf("<strong>PART 1 </strong>\n");
% Ch. 1, Problem 1
d = dms2degrees([0,1,0]); % d [deg] = 1 [min] = 1852 [m]
d1 = d/1852*0.01; % linearize d [deg] to 0.01 [m]
d2 = 1*60/1852*0.01; % linearize d [sec] to 0.01 [m]
fprintf("1a) In decimal degrees, the decimal precision is 10^%d. \n", floor(log10(d1)));
fprintf("1b) In DMS degrees, the decimal precision is 10^%d. \n", floor(log10(d2)));

% Ch. 1, Problem 2
llh = [45, -120, 10e3]; % [deg, deg, m]
v = 885e3/3600; % [m/s]
az = 45; % [deg]
t = 8*3600; % [s]
R = 6371e3; % [m] - earth radius
% w = 15; % [deg/hr] - earth rotation

xyz = lla2enu(llh, llh, "flat");
dx = [sind(az), cosd(az), 0] .* v .* t; % [m]
xyz2 = xyz + dx;
llh2 = enu2lla(xyz2, llh, "flat");

vn = v * sind(az);
ve = v * cosd(az);

for i = 1:t
    if i == 1
        lat(i) = deg2rad(llh(1));
        lon(i) = deg2rad(llh(2));
        hgt(i) = llh(3);
    else
        % groves 2.111
        dt = 1;
        dLat = vn / (R + hgt(i-1));
        dLon = ve / ((R + hgt(i-1)) * cos(lat(i-1)));

        lat(i) = lat(i-1) + dt*dLat;
        lon(i) = lon(i-1) + dt*dLon;
        hgt(i) = hgt(i-1) + dt*0;
    end
end

f = uifigure(HandleVisibility="on");
g = geoglobe(f);
geoplot3(g, rad2deg(lat), rad2deg(lon), hgt, 'c', LineWidth=2);

fprintf("2) New LLH = [%f, %f, %f] \n", llh2(1), llh2(2), llh2(3)/1000);
fprintf("2) New LLH = [%f, %f, %f] \n", rad2deg(lat(end)), rad2deg(lon(end)), hgt(end)/1000);

% Ch. 1, Problem 3
dt = 0.1; % [s]
v = 360e3/3600; % [m/s]
f = 100e6; % [Hz]
L = 3e8 / f;
df = [-33.1679; -33.1711; -33.1743]; % [Hz]
r_dot = -L .* df;
fprintf("3a) The range rate for each doppler shift is [%f, %f, %f] \n", r_dot);

% syms x0 y0 V
% x1 = x0 + 10;
% x2 = x0 + 20;
% 
% f = [x0 / sqrt(x0^2 + y0^2);...
%      x1 / sqrt(x1^2 + y0^2); ...
%      x2 / sqrt(x2^2 + y0^2)] .* v; %.* v;
% J = simplify(jacobian(f, [x0,y0]));
% H_func = @(x,y) double(subs(J, [x0, y0],[x, y]));
% y_func = @(x,y) r_dot - double(subs(f, [x0, y0], [x, y]));

% x = fsolve(@myfunc, [1,1])
% fprintf("%d, %f \n", x)

% est = [1;1];
% err = 10;
% l = 10;
% while err > 1e-12
%     H = H_func(est(1), est(2));
%     y = y_func(est(1), est(2));
% 
%     del = inv(H'*H)*H' * y;
%     est = est + del;
% %     t_err = dot(del, del);
%     err = norm(del);
% % 
% %     if t_err < err
% %         l = l/10;
% %         est = t_est;
% %         err = t_err;
% %     else
% %         l = l*10;
% %     end
% end
% est

% Ch. 1, Problem 4
syms x b;
p1 = 550 == sqrt((0-x)^2) - b;
p2 = 500 == sqrt((1000-x)^2) - b;
[b1, x1] = solve(p1,p2);
fprintf("4a) Given p1 = 550 and p2 = 500, x = %f and b = %f \n", x1, b1);
syms x b;
p1 = 400 == sqrt(x^2) - b;
p2 = 1400 == sqrt(1000-x)^2 - b;
[b2, x2] = solve(p1,p2);
fprintf("4b) Given p1 = 400 and p2 = 1400, x = %f and b = %f \n", x2, b2);

fprintf("\n");

%% PART 2
fprintf("<strong>PART 2 </strong>\n");
fprintf("a bunch of plots! \n")
fprintf("\n");

seq1 = 2*ceil(rand(100,1)-0.5)-1;
seq1_auto = correlation(seq1, seq1);
seq1_xcorr = xcorr(seq1, seq1, 'normalized');
seq1_1000 = 2*ceil(rand(1000,1)-0.5)-1;
seq1_1000_auto = correlation(seq1_1000, seq1_1000);
% seq1_1000_xcorr = xcorr(seq1_1000, seq1_1000);

seq2 = 2*ceil(0.1*randn(100,1))-1;
seq2_auto = correlation(seq2, seq2);
% seq2_xcorr = xcorr(seq2,seq2);
seq2_1000 = 2*ceil(0.1*randn(1000,1))-1;
seq2_1000_auto = correlation(seq2_1000, seq2_1000);
% seq2_1000_xcorr = xcorr(seq2, seq2);

cross_corr = correlation(seq1, seq2);
cross_corr_1000 = correlation(seq1_1000, seq2_1000);
% cross_corr = xcorr(seq1, seq2);
% cross_corr_1000 = xcorr(seq1_1000, seq2_1000);

f = figure(Name="PART 2");
tabs = uitabgroup(f);
tab1 = uitab(tabs, Title="seq1 - rand");
tab2 = uitab(tabs, Title="seq2 - randn");
tab3 = uitab(tabs, Title="Cross Correlation");
tab4 = uitab(tabs, Title="seq1_1000 - rand");
tab5 = uitab(tabs, Title="seq2_1000 - randn");
tab6 = uitab(tabs, Title="Cross Correlation 1000");

axes(Parent=tab1);
subplot(3,1,1);
histogram(seq1, FaceColor="red");
grid("on");
title("Histogram");
subplot(3,1,2);
% periodogram(seq1);
% pwelch(seq1,window_filter);
plot(abs(fft(seq1)), Color="red");
grid("on");
title("Power Spectral Density");
subplot(3,1,3);
hold on
stem(-99:99, seq1_auto, Color="red");
stem(-99:99, seq1_xcorr, Color="cyan");
hold off
grid("on");
title("Auto Correlation");

axes(Parent=tab2);
subplot(3,1,1);
histogram(seq2, FaceColor="blue");
grid("on");
title("Histogram");
subplot(3,1,2);
% periodogram(seq2);
% pwelch(seq2,window_filter);
plot(abs(fft(seq2)), Color="blue");
grid("on");
title("Power Spectral Density");
subplot(3,1,3);
stem(-99:99, seq2_auto, Color="blue");
% stem(seq2_xcorr, Color="red");
title("Auto Correlation");

axes(Parent=tab3);
subplot(3,1,2);
stem(-99:99, cross_corr, 'g');
% stem(cross_corr, Color="green");
grid("on");
title("Correlation between the Sequences");

axes(Parent=tab4);
subplot(3,1,1);
histogram(seq1_1000, FaceColor="red");
grid("on");
title("Histogram");
subplot(3,1,2);
% periodogram(seq1);
% pwelch(seq1,window_filter);
plot(abs(fft(seq1_1000)), Color="red");
grid("on");
title("Power Spectral Density");
subplot(3,1,3);
stem(-999:999, seq1_1000_auto, Color="red");
% stem(seq1_1000_xcorr, Color="red");
grid("on");
title("Auto Correlation");

axes(Parent=tab5);
subplot(3,1,1);
histogram(seq2_1000, FaceColor="blue");
grid("on");
title("Histogram");
subplot(3,1,2);
% periodogram(seq2);
% pwelch(seq2,window_filter);
plot(abs(fft(seq2_1000)), Color="blue");
grid("on");
title("Power Spectral Density");
subplot(3,1,3);
stem(-999:999, seq2_1000_auto, Color="blue");
% stem(seq2_1000_xcorr, Color="red");
title("Auto Correlation");

axes(Parent=tab6);
subplot(3,1,2);
stem(-999:999, cross_corr_1000, 'g');
% stem(cross_corr_1000, Color="green");
grid("on");
title("Correlation between the Sequences");

%% PART 3
fprintf("<strong>PART 3 </strong>\n");

A = 3 + 3*randn(1000,1);
B = 5 + 5*randn(1000,1);
C = A + B;
DATA = [A, B, C];

fprintf("Mean: A = %f \n      B = %f \n      C = %f \n", mean(A), mean(B), mean(C));
fprintf("StdDev: A = %f \n        B = %f \n        C = %f \n", std(A), std(B), std(C));
fprintf("Variance: A = %f \n          B = %f \n          C = %f \n", var(A), var(B), var(C));
fprintf("\n");

fprintf("Mean of DATA = [A B C] = %f \n", mean(DATA, "all"));
fprintf("Covariance of DATA = \n");
disp(cov(DATA));

fprintf("\n");

%% PART 4
fprintf("<strong>PART 4 </strong>\n");

syms x y a b x0 y0
r = sqrt( (x-a)^2 + (y-b)^2 );
pretty(taylor(r, [x,y], ExpansionPoint=[x0,y0], Order=2));

fprintf("\n");


%% FUNCTIONS
function [R, shifts] = correlation(seq1, seq2)
    if size(seq1) ~= size(seq2)
        fprintf("ERROR: 'correlation' -> sequneces not the same size!");
        return;
    end

    N = length(seq1);
    shifts = -(N-1):(N-1);
    R = zeros([length(shifts),1]);

    for k = 1:length(shifts)
        seq2_ = circshift(seq2, shifts(k));
        R(k) = 1/N * sum(seq1.*seq2_);
    end
end

function F = myfunc(x)
x0 = x(1);
y0 = x(2);
x1 = x0 + 10;
x2 = x0 + 20;
v = 100;
F = [x0 / sqrt(x0^2 + y0^2);...
     x1 / sqrt(x1^2 + y0^2); ...
     x2 / sqrt(x2^2 + y0^2)] .* v; %.* v;
end