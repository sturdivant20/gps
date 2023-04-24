%% GPS HW3 - Problem 4 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 4</strong>\n");

% bevly generated signal
signal = generate_signal(3);
code = [1, -1, -1, -1, 1, -1, 1, 1];
codeLen = length(code);
nSamples = 80;
prn = repelem(code, 10);

% timing
t_int = 1;
f_s = 80/t_int;
t_s = 1/f_s;

% correlator spacing
s = 5; % number of samples

% % controller design
% BW = 1;                 % controller bandwidth (1-5 Hz)
% K = 1;                  % loop gain
% z = 2*sind(45);         % damping ratio
% w = BW*4*z / (1+z^2);   % natural frequency
% 
% % PI
% Ki = w^2;
% Kp = z*w;
Kp = 5;

fig = figure(Units='normalized', Position=[3.0, 0.4, 1.2, 0.4]);
loop = 1:80:length(signal);
L = length(loop);
disc = zeros(1,L);
tau = zeros(1,L);

beg = 1;
fin = nSamples;

% upsampling index
cSamp = 10; % samples per chip
j = 0;
idx = ceil(1/nSamples : 8/nSamples : 8-(1/nSamples));

for i = 2:L

    % signal chunk
    sig = signal(beg:fin);
    
    % autocorrelatation
    E = sig * circshift(prn, tau(i-1)+5)';
    P = sig * circshift(prn, tau(i-1))';
    L = sig * circshift(prn, tau(i-1)-5)';

    % discriminator
    disc(i) = 0.5 * (E-L) / (E+L);

    % propagate
    tau(i) = floor(tau(i-1) + Kp*disc(i));

    beg = fin + 1;
    fin = fin + 80;

%     % gif
%     if i < 3
%         if exist("./media/part4.gif")
%             delete("./media/part4.gif")
%         end
%     end
%     tl = tiledlayout(2,1, TileSpacing="tight");
%     ax = axes(Parent=tl);
%     hold on;
%     stairs(sig, 'b', LineWidth=2)
%     stairs(prn, 'r--', LineWidth=2);
%     legend("Signal", "Shifted PRN");
%     title("Shifted Signal");
%     ylim([-1.1,1.1])
%     hold off;
%     nexttile;
%     stairs(disc, 'b', LineWidth=2);
%     title("Discriminator Error")
%     ylim([-0.6,0.6])
%     exportgraphics(fig, './media/part4.gif', Append=true)

end

figure;
tl = tiledlayout(3,1, TileSpacing="tight");
ax = axes(Parent=tl);
hold on;
stairs(sig, 'b', LineWidth=2)
stairs(circshift(prn, tau(i-1)), 'r--', LineWidth=2);
legend("Signal", "Shifted PRN", Location="southoutside", Orientation="horizontal");
title("Shifted Signal");
ylim([-1.1,1.1])
hold off;
nexttile;
stairs(disc, 'b', LineWidth=2);
title("Discriminator Error")
ylim([-0.6,0.6])
nexttile;
stairs(tau, 'b', LineWidth=2)
title("PRN Ofset (\tau)")
