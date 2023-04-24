%% GPS HW3 - Problem 5 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 5</strong>\n");

% f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
% tbs = uitabgroup(Parent=f);
% tab(1) = uitab(Parent=tbs, Title="Autocorrelation");
% tab(2) = uitab(Parent=tbs, Title="noisy Auocorrelation");

load("ca_codes.mat")
PRN4 = repelem(ca_code(4,:),16);
PRN7 = repelem(ca_code(7,:),16);

shifts = -5*16:5*16;
R = correlation(PRN4, PRN4, 'standard', shifts);
R_noise = correlation(PRN4+0.2*randn(size(PRN4)), PRN4, 'standard', shifts);

figure
tiledlayout(2,1, TileSpacing="tight");
nexttile;
stem(shifts, R, 'r', LineWidth=2)
title("Perfect AutoCorrelation")
nexttile;
stem(shifts, R_noise, 'r', LineWidth=2)
title("Noisy AutoCorrelation")
xlabel("Bit Shifts")