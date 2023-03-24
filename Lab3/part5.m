%% GPS LAB 3 - PART 5 | Daniel Sturdivant & Andrew Weir
clc; clear; close all;
fprintf("<strong>PART 4\n</strong>");

load("RCVR_S1.mat");
load("RCVR_D1.mat");

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

clearvars -except S1 D1 
L = length(D1);


%% Part A

% ephem = utils.readMultiGNSSRinex("part5Data.rnx", 2253);
[ephem, ~, ~] = utils.ReadNav("part5Data.rnx");
% everything is rinex 2?
