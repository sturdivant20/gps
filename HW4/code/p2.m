%% GPS HW3 - Problem 2 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 2</strong>\n");

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
f = figure(Units='normalized', Position=[3.0, 0.4, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Parent=tbs, Title="A) PLL");
tab(2) = uitab(Parent=tbs, Title="A) PLL_IQ");
tab(3) = uitab(Parent=tbs, Title="A) PLL_wave");
tab(4) = uitab(Parent=tbs, Title="B) PLL");
tab(5) = uitab(Parent=tbs, Title="B) PLL_IQ");
tab(6) = uitab(Parent=tbs, Title="B) PLL_wave");
tab(7) = uitab(Parent=tbs, Title="C) PLL");
tab(8) = uitab(Parent=tbs, Title="C) PLL_IQ");
tab(9) = uitab(Parent=tbs, Title="C) PLL_wave");
tab(10) = uitab(Parent=tbs, Title="D) PLL");
tab(11) = uitab(Parent=tbs, Title="D) PLL_IQ");
tab(12) = uitab(Parent=tbs, Title="D) PLL_wave");

% bevly encrypted signal
signal = generate_signal(1);
L = length(signal);


%% Part A
fprintf("\n<strong>(a)</strong>\n");

% controller design
BW = 1;                 % controller bandwidth
K = 4;                  % loop gain
z = sind(45);           % damping ratio
w = BW*4*z / (1+z^2);   % natural frequency

% PI
Ki = w^2;
Kp = z*w;

% time setup
t_s = 1e-6;             % sampling time
t_int = 0.01;           % integration time
t = 0:t_s:t_int-t_s;        % simulated signal time
numSamp = t_int / t_s;  % samples per integrastion period
loop = 1:numSamp:L;

% initlize values
disc = zeros(size(loop));   % error
theta = zeros(size(loop));  % phase angle
f = zeros(size(loop));      % frequency
I = zeros(size(loop));      % Inphase
Q = zeros(size(loop));      % Quadrature
f(1) = 100;

beg = 1;
fin = numSamp;

axes(parent=tab(3));
for i = 2:length(loop)

    % data signal
    sig = signal(beg:fin);

    % Inphase and Quadrature
    ssig = sin(2*pi*f(i-1).*t' + theta(i-1));
    csig = cos(2*pi*f(i-1).*t' + theta(i-1));
    I(i) = sig * ssig;
    Q(i) = sig * csig;

    % error
    disc(i) = atan(Q(i) / I(i));

    % propagate
    f(i) = f(i-1) + K * (Ki*t_int*disc(i) + Kp*(disc(i) - disc(i-1)));
    theta(i) = rem(theta(i-1) + 2*pi*f(i-1)*t_int, 2*pi);

    beg = beg + numSamp;
    fin = fin + numSamp;


%     % gif
%     if i < 2
%         if exist("./media/part2_a.gif")
%             delete("./media/part2_a.gif")
%         end
%     end
%     plot(sig, 'r')
%     hold on;
%     plot(ssig, 'b')
%     legend("Data", "Sin")
%     hold off;
%     exportgraphics(tab(3), './media/part2_a.gif', Append=true)

end

tl = tiledlayout(3,1, Parent=tab(1), TileSpacing="tight");
ax =axes(Parent=tl);
plot(disc)
title("Error")
xlim([0,101])

nexttile;
plot(theta)
title("Phase")
xlim([0,101])

nexttile;
hold on;
plot(f);
% plot(eHat+100);
title("Frequency")
% legend("Frequency", "Controller Error+100");
xlim([0,101])

ax = axes(Parent=tab(2));
hold on
plot(I)
plot(Q)
legend("Inphase", "Quadurature");
xlim([0,101])

ax = axes(Parent=tab(3));
hold on
plot(sig, 'r')
plot(ssig, 'b')
title("Signal Comparison at Final Iteration")
legend("Actual Signal", "Simulated Signal");


%% Part B
fprintf("\n<strong>(b)</strong>\n");

% controller design
BW = 2;                 % controller bandwidth
K = 4;                  % loop gain
z = sind(45);           % damping ratio
w = BW*4*z / (1+z^2);   % natural frequency

% PI
Ki = w^2;
Kp = z*w;

% time setup
t_s = 1e-6;             % sampling time
t_int = 0.01;           % integration time
t = 0:t_s:t_int-t_s;        % simulated signal time
numSamp = t_int / t_s;  % samples per integrastion period
loop = 1:numSamp:L;

% initlize values
disc = zeros(size(loop));    % error
theta = zeros(size(loop));  % phase angle
f = zeros(size(loop));      % frequency
I = zeros(size(loop));      % Inphase
Q = zeros(size(loop));      % Quadrature
f(1) = 100;

beg = 1;
fin = numSamp;

axes(parent=tab(6));
for i = 2:length(loop)

    % data signal
    sig = signal(beg:fin);

    % Inphase and Quadrature
    ssig = sin(2*pi*f(i-1).*t' + theta(i-1));
    csig = cos(2*pi*f(i-1).*t' + theta(i-1));
    I(i) = sig * ssig;
    Q(i) = sig * csig;

    % error
    disc(i) = atan2(Q(i) , I(i));

    % propagate
    f(i) = f(i-1) + K * (Ki*t_int*disc(i) + Kp*(disc(i) - disc(i-1)));
    theta(i) = rem(theta(i-1) + 2*pi*f(i-1)*t_int, 2*pi);

    beg = beg + numSamp;
    fin = fin + numSamp;


%     % gif
%     if i < 2
%         if exist("./media/part2_b.gif")
%             delete("./media/part2_b.gif")
%         end
%     end
%     plot(sig, 'r')
%     hold on;
%     plot(ssig, 'b')
%     legend("Data", "Sin")
%     hold off;
%     exportgraphics(tab(6), './media/part2_b.gif', Append=true)

end

tl = tiledlayout(3,1, Parent=tab(4), TileSpacing="tight");
ax =axes(Parent=tl);
plot(disc)
title("Error")
xlim([0,101])

nexttile;
plot(theta)
title("Phase")
xlim([0,101])

nexttile;
hold on;
plot(f);
% plot(eHat+100)
title("Frequency")
% legend("Frequency", "Controller Error+100");
xlim([0,101])

ax = axes(Parent=tab(5));
hold on
plot(I)
plot(Q)
legend("Inphase", "Quadurature");
xlim([0,101])

ax = axes(Parent=tab(6));
hold on
plot(sig, 'r')
plot(ssig, 'b')
title("Signal Comparison at Final Iteration")
legend("Actual Signal", "Simulated Signal");


%% Part C
fprintf("\n<strong>(c)</strong>\n");

% controller design
BW = 2;                 % controller bandwidth
K = 4;                  % loop gain
z = sind(45);           % damping ratio
w = BW*4*z / (1+z^2);   % natural frequency

% PI
Ki = w^2;
Kp = z*w;

% time setup
t_s = 1e-6;             % sampling time
t_int = 0.01;           % integration time
t = 0:t_s:t_int-t_s;        % simulated signal time
numSamp = t_int / t_s;  % samples per integrastion period
loop = 1:numSamp:L;

% initlize values
disc = zeros(size(loop));    % error
theta = zeros(size(loop));  % phase angle
f = zeros(size(loop));      % frequency
I = zeros(size(loop));      % Inphase
Q = zeros(size(loop));      % Quadrature
f(1) = 110;

beg = 1;
fin = numSamp;

axes(parent=tab(9));
for i = 2:length(loop)

    % data signal
    sig = signal(beg:fin);

    % Inphase and Quadrature
    ssig = sin(2*pi*f(i-1).*t' + theta(i-1));
    csig = cos(2*pi*f(i-1).*t' + theta(i-1));
    I(i) = sig * ssig;
    Q(i) = sig * csig;

    % error
    disc(i) = atan(Q(i) / I(i));

    % propagate
    f(i) = f(i-1) + K * (Ki*t_int*disc(i) + Kp*(disc(i) - disc(i-1)));
    theta(i) = rem(theta(i-1) + 2*pi*f(i-1)*t_int, 2*pi);

    beg = beg + numSamp;
    fin = fin + numSamp;


%     % gif
%     if i < 2
%         if exist("./media/part2_c.gif")
%             delete("./media/part2_c.gif")
%         end
%     end
%     plot(sig, 'r')
%     hold on;
%     plot(ssig, 'b')
%     legend("Data", "Sin")
%     hold off;
%     exportgraphics(tab(9), './media/part2_c.gif', Append=true)

end

tl = tiledlayout(3,1, Parent=tab(7), TileSpacing="tight");
ax =axes(Parent=tl);
plot(disc)
title("Error")
xlim([0,101])

nexttile;
plot(theta)
title("Phase")
xlim([0,101])

nexttile;
hold on;
plot(f);
% plot(eHat+110)
title("Frequency")
% legend("Frequency", "Controller Error+110");
xlim([0,101])

ax = axes(Parent=tab(8));
hold on
plot(I)
plot(Q)
legend("Inphase", "Quadurature");
xlim([0,101])

ax = axes(Parent=tab(9));
hold on
plot(sig, 'r')
plot(ssig, 'b')
title("Signal Comparison at Final Iteration")
legend("Actual Signal", "Simulated Signal");


%% Part D
fprintf("\n<strong>(d)</strong>\n");

% controller design
BW = 5;                 % controller bandwidth
K = 4;                  % loop gain
z = sind(45);           % damping ratio
w = BW*4*z / (1+z^2);   % natural frequency

% PI
Ki = w^2;
Kp = z*w;

% time setup
t_s = 1e-6;             % sampling time
t_int = 0.01;           % integration time
t = 0:t_s:t_int-t_s;        % simulated signal time
numSamp = t_int / t_s;  % samples per integrastion period
loop = 1:numSamp:L;

% initlize values
disc = zeros(size(loop));    % error
theta = zeros(size(loop));  % phase angle
f = zeros(size(loop));      % frequency
I = zeros(size(loop));      % Inphase
Q = zeros(size(loop));      % Quadrature
f(1) = 10;

beg = 1;
fin = numSamp;

axes(parent=tab(12));
for i = 2:length(loop)

    % data signal
    sig = signal(beg:fin);

    % Inphase and Quadrature
    ssig = sin(2*pi*f(i-1).*t' + theta(i-1));
    csig = cos(2*pi*f(i-1).*t' + theta(i-1));
    I(i) = sig * ssig;
    Q(i) = sig * csig;

    % error
    disc(i) = atan(Q(i) / I(i));

    % propagate
    f(i) = abs(f(i-1) + K * (Ki*t_int*disc(i) + Kp*(disc(i) - disc(i-1))));
    theta(i) = rem(theta(i-1) + 2*pi*f(i-1)*t_int, 2*pi);

    beg = beg + numSamp;
    fin = fin + numSamp;


%     % gif
%     if i < 2
%         if exist("./media/part2_d.gif")
%             delete("./media/part2_d.gif")
%         end
%     end
%     plot(sig, 'r')
%     hold on;
%     plot(ssig, 'b')
%     legend("Data", "Sin")
%     hold off;
%     exportgraphics(tab(9), './media/part2_d.gif', Append=true)

end

tl = tiledlayout(3,1, Parent=tab(10), TileSpacing="tight");
ax =axes(Parent=tl);
plot(disc)
title("Error")
xlim([0,101])

nexttile;
plot(theta)
title("Phase")
xlim([0,101])

nexttile;
hold on;
plot(f);
% plot(eHat+10)
title("Frequency")
% legend("Frequency", "Controller Error+10");
xlim([0,101])

ax = axes(Parent=tab(11));
hold on
plot(I)
plot(Q)
legend("Inphase", "Quadurature");
xlim([0,101])

ax = axes(Parent=tab(12));
hold on
plot(sig, 'r')
plot(ssig, 'b')
title("Signal Comparison at Final Iteration")
legend("Actual Signal", "Simulated Signal");

%%

exportgraphics(tab(1), './media/p2_a.png')
exportgraphics(tab(3), './media/p2_a_wave.png')

exportgraphics(tab(4), './media/p2_b.png')
exportgraphics(tab(6), './media/p2_b_wave.png')

exportgraphics(tab(7), './media/p2_c.png')
exportgraphics(tab(9), './media/p2_c_wave.png')

exportgraphics(tab(10), './media/p2_d.png')
exportgraphics(tab(12), './media/p2_d_wave.png')

