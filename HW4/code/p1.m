%% GPS HW3 - Problem 1 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 1</strong>\n");

f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
% f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Parent=tbs, Title="Unit Input");
tab(2) = uitab(Parent=tbs, Title="Ramp Input");


%% Part A
fprintf("\n<strong>(a)</strong>\n");

ts = 1;         % settling time
z = sind(45);    % damping
w = 4.6/(z*ts); % natural frequency

% % PID Controller
% Kd = 1;
% Ki = w^2 * Kd;
% Kp = 2*z*w * Kd - 1;
% 
% K = tf([Kd, Kp, Ki], [1, 0]);
% G = tf(1, [1,0]);

% PI Controller
Ki = w^2;
Kp = 2*z*w;

K = tf([Kp, Ki], [1, 0]);
G = tf(1, [1,0]);

% feedback system
CL = minreal((K*G)/(1 + K*G));
s = eig(CL);
[Gm,Pm] = margin(CL);

% ss error to unit input
time = 0:0.01:3;
r = ones(size(time));
y = lsim(CL, r, time);

ax = axes(Parent=tab(1));
hold on;
plot(time, r, Color=[0.5 0.5 0.5]);
plot(time, y, 'b');
title("Linear Simulation Results")
xlabel("Time [seconds]")
ylabel("Amplitude")


%% Part B
fprintf("\n<strong>(b)</strong>\n");

% ss error to ramp input
time = 0:0.01:3;
r = 1:0.01:4;
y = lsim(CL, r, time);

ax = axes(Parent=tab(2));
hold on;
plot(time, r, Color=[0.5 0.5 0.5]);
plot(time, y, 'b');
title("Linear Simulation Results")
xlabel("Time [seconds]")
ylabel("Amplitude")


%% Part C
fprintf("\n<strong>(c)</strong>\n");

%%
figure 
margin(CL)

