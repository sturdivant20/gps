%% Daniel Sturdivant || MECH 6970 GPS HW3
clc; clear; close all;

% f = figure('units','normalized','position',[0.1 0.1 0.3 0.7]);
% f = figure('units','normalized','position',[0.1 0.1 0.7 0.8]);
f = figure('units','normalized','position',[0.1 0.1 0.9 1.0]);
tabs = uitabgroup(f);
tab(1) = uitab(Title='1b)');
tab(2) = uitab(Title='1c)');
tab(3) = uitab(Title='4)');
tab(4) = uitab(Title='5)');
tab(5) = uitab(Title='6a)');
tab(6) = uitab(Title='6b)');
tab(7) = uitab(Title='7a)');
tab(8) = uitab(Title='7b)');
tab(9) = uitab(Title='7c)');
tab(10) = uitab(Title='8)');

%% PROBLEM 1
fprintf("\n<strong>PROBLEM 1</strong> \n");

% PART A
fprintf("(a)\n");

std_a = 1;
std_b = 1;

% what standard deviation in y should be
std_y = sqrt(9*std_a^2 + 16*std_b^2);

% 1000 run monte carlo to confirm
y1 = zeros(1,1000);
y2 = zeros(1,1000);
for i = 1:1000
    a = 1 + std_a * randn;
    b = 1 + std_b * randn;

    y1(i) = 3*a + 4*b;
    y2(i) = 3*a - 4*b;
end

std_y1 = std(y1);
std_y2 = std(y2);

fprintf("Expected sigma_y = %f\n", std_y);
fprintf("MC Sum   sigma_y = %f\n", std_y1);
fprintf("MC Diff. sigma_y = %f\n", std_y2);

% PART B
fprintf("(b)\n");

time = 0:10*60;
dt = 1;
noise = [0.01, 0.1];

w = zeros(1000, length(time), 2);
sig = zeros(length(time), 2);
sig_mc = zeros(length(time), 2);
mu = zeros(length(time), 2);
mu_mc = zeros(length(time), 2);

% 2 standard deviations
for n = 1:2
    % 1000 MC runs
    for m = 1:1000
        % 10 minute run
        for t = 1:length(time)
            if t == 1
                w(m,t,n) = noise(n)*randn;
                sig(t,n) = noise(n) * dt * sqrt(t);
            else
                w(m,t,n) = w(m,t-1,n) + noise(n)*randn;
                sig(t,n) = noise(n) * dt * sqrt(t);
            end
        end
    end
end
sig_mc(:,1) = std(w(:,:,1));
sig_mc(:,2) = std(w(:,:,2));
mu_mc(:,1) = mean(w(:,:,1));
mu_mc(:,2) = mean(w(:,:,2));

tl = tiledlayout(2,1, TileSpacing="tight", Parent=tab(1));
title(tl, "\textbf{1000 Run Monte Carlo of Integrated White Noise}", Interpreter="latex", FontSize=18);
xlabel(tl, "Time [s]", FontSize=16);
ylabel(tl, "White Noise", FontSize=16);
ax(1) = axes(Parent=tl);
ax(1).XAxis.Visible = "off";
hold on;
plot(time, w(:,:,1), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,1), Color="#EC142F", LineWidth=2.5, DisplayName="\pm1\sigma=0.01");
plot(time, -sig(:,1), Color="#EC142F", LineWidth=2.5, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="Northwest");
nexttile;
hold on;
plot(time, w(:,:,2), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,2), Color="#EC142F", LineWidth=2.5, DisplayName="\pm1\sigma=0.1");
plot(time, -sig(:,2), Color="#EC142F", LineWidth=2.5, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="Northwest");
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% PART C
fprintf("(c)\n");

time = 0:10*60;
dt = 1;
noise = [0.01, 0.1];
time_const = [1, 100];

x = zeros(1000, length(time), 4);
sig = zeros(length(time), 4);
sig_mc = zeros(length(time), 4);
mu = zeros(length(time), 4);
mu_mc = zeros(length(time), 4);

% 2 standard deviations
for i = 1:2
    n = noise(i);

    % 2 time constants
    for s = 1:2
        tau = time_const(s);
        if i == 1
            ii = i + s - 1;
        else
            ii = i + s;
        end

        % 1000 MC runs
        for m = 1:1000

            % 10 minute run
            for t = 1:length(time)
                w = n*randn;
%                 x(m,t,ii) = w*tau*(1 + exp(-time(t)/tau));
                if t == 1
                    x(m,t,ii) = w*dt;
                else
                    xDot = -x(m,t-1,ii)/tau + w;
                    x(m,t,ii) = x(m,t-1,ii) + xDot*dt;
                end
                A = 1 - (dt / tau);
                sig(t,ii) = n * dt * sqrt((A^(2*time(t)) - 1) / (A^2 - 1));
            end
        end
    end
end

tl = tiledlayout(4,1, TileSpacing="tight", Parent=tab(2));
title(tl, "\textbf{1000 Run Monte Carlo of Markov Process}", Interpreter="latex", FontSize=18);
xlabel(tl, "Time [s]", FontSize=16);
ylabel(tl, "Markov Noise", FontSize=16)
a1 = axes(Parent=tl);
a1.XAxis.Visible = "off";
hold on;
plot(time, x(:,:,1), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,1), Color="#EC142F", LineWidth=3, DisplayName="\pm1\sigma=0.01");
plot(time, -sig(:,1), Color="#EC142F", LineWidth=3, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="North", Orientation="Horizontal", FontSize=14);
title("\tau=1", FontSize=16);
a3 = nexttile;
a3.XAxis.Visible = "off";
hold on;
plot(time, x(:,:,3), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,3), Color="#EC142F", LineWidth=3, DisplayName="\pm1\sigma=0.1");
plot(time, -sig(:,3), Color="#EC142F", LineWidth=3, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="North", Orientation="Horizontal", FontSize=14);
a2 = nexttile;
a2.XAxis.Visible = "off";
hold on;
plot(time, x(:,:,2), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,2), Color="#EC142F", LineWidth=3, DisplayName="\pm1\sigma=0.01");
plot(time, -sig(:,2), Color="#EC142F", LineWidth=3, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="North", Orientation="Horizontal", FontSize=14);
title("\tau=100", FontSize=16);
a4 = nexttile;
hold on;
plot(time, x(:,:,4), Color="#36454F", HandleVisibility="off");
plot(time, sig(:,4), Color="#EC142F", LineWidth=3, DisplayName="\pm1\sigma=0.1");
plot(time, -sig(:,4), Color="#EC142F", LineWidth=3, HandleVisibility="off");
plot(time, zeros(size(time)), Color="#EDB120", LineWidth=2, DisplayName="\mu=0");
grid on;
legend(Location="North", Orientation="Horizontal", FontSize=14);
linkaxes([a1, a2, a3, a4], 'x');


%% PROBLEM 2
fprintf("\n<strong>PROBLEM 2</strong> \n");

f1 = 1575.42e6;
f_L2 = 1227.60e6;
f_L5 = 1176.45e6;
syms rho_L1 rho_L2 rho_L5 real

psr_IF_L1L2 = (f1^2 / (f1^2 - f_L2^2))*rho_L1 + (f_L2^2 / (f1^2 - f_L2^2))*rho_L2;
psr_IF_L1L5 = (f1^2 / (f1^2 - f_L5^2))*rho_L1 + (f_L5^2 / (f1^2 - f_L5^2))*rho_L5;
psr_IF_L2L5 = (f_L2^2 / (f_L2^2 - f_L5^2))*rho_L2 + (f_L5^2 / (f_L2^2 - f_L5^2))*rho_L5;

fprintf("psr_IF_L1L2 = %s\n", char(psr_IF_L1L2));
fprintf("psr_IF_L1L5 = %s\n", char(psr_IF_L1L5));
fprintf("psr_IF_L2L5 = %s\n", char(psr_IF_L2L5));

sig_L1L2 = sqrt((5929/2329)^2 + (3600/2329)^2);
sig_L1L5 = sqrt((23716/10491)^2 + (13225/10491)^2);
sig_L2L5 = sqrt((576/47)^2 + (576/47)^2);

fprintf("sig_L1L2 = %f\n", sig_L1L2);
fprintf("sig_L1L5 = %f\n", sig_L1L5);
fprintf("sig_L2L5 = %f\n", sig_L2L5);


%% PROBLEM 3
fprintf("\n<strong>PROBLEM 3</strong> \n");


%% PROBLEM 4
fprintf("\n<strong>PROBLEM 4</strong> \n");

x_sv = [  0, 300; ...
        100, 400; ...
        700, 400; ...
        800, 300];
x_r = [400, 0];
x_u = [401, 0];

R_r = sqrt(sum((x_sv-x_r).^2, 2));
R_u = sqrt(sum((x_sv-x_u).^2, 2));

% PART A 
x1 = [0;0];
e1 = 1;
x2 = [0;0];
e2 = 1;
while (e1 > 1e-6)
    % 2 sv solution
    u1 = x_sv(1:2,:) - x1';
    R1 = sqrt(sum(u1.^2,2));
    uv1 = u1./R1;

    H1 = -uv1;
    y1 = R_u(1:2,:) - R1;

    dx1 = pinv(H1)*y1;
    x1 = x1 + dx1;
    e1 = norm(dx1);
end
while (e2 > 1e-6)
    % 4 sv solution
    u2 = x_sv - x2';
    R2 = sqrt(sum(u2.^2,2));
    uv2 = u2./R2;

    H2 = -uv2;
    y2 = R_u - R2;

    dx2 = pinv(H2)*y2;
    x2 = x2 + dx2;
    e2 = norm(dx2);
end
DOP1 = inv(H1'*H1);
DOP2 = inv(H2'*H2);

fprintf("(a)\nPDOP1 = %f\n", sqrt(trace(DOP1)));
fprintf("PDOP2 = %f\n", sqrt(trace(DOP2)));

% PART B
x2 = [0;0;0];
e2 = 1;
while (e2 > 1e-6)
    u2 = x_sv - x2(1:2)';
    R2 = sqrt(sum(u2.^2,2));
    uv2 = u2./R2;

    H2 = [-uv2, ones(4,1)];
    y2 = R_u - (R2 + x2(3));
    dx2 = pinv(H2)*y2;

    x2 = x2 + dx2;
    e2 = norm(dx2);
end
DOP2 = inv(H2'*H2);

fprintf("(b)\nPDOP2 = %f\n", sqrt(trace(DOP2(1:2,1:2))));
fprintf("GDOP2 = %f\n", sqrt(trace(DOP2)));

% PART C
H3 = [(x_sv-x_r)./R_r, ones(4,1)];
y3 = R_r - R_u;
dx3 = pinv(H3)*y3;
x3 = [x_r,0]' + dx3;
DOP3 = inv(H3'*H3);

fprintf("(c)\nPDOP = %f\n", sqrt(DOP3(1,1)+DOP3(2,2)));
fprintf("GDOP = %f\n", sqrt(trace(DOP3)));

% PART D
u = (x_sv-x_r)./R_r;
dR = R_r - R_u;
H4 = u(2:4,:) - u(1,:);
y4 = dR(2:4) - dR(1);
dx4 = pinv(H4)*y4;
x4 = x_r' + dx4;
DOP4 = inv(H4'*H4);

fprintf("(d)\nPDOP = %f\n", sqrt(trace(DOP4)));


%% PROBLEM 5
fprintf("\n<strong>PROBLEM 5</strong> \n");

% 32.587, -85.4939, 183
gpsTime = 523323; % i=121
svInUse = [3; 6; 11; 14; 17; 19; 24; 30];
psr1 = [23223192.3246802; ...
        20670605.644634; ...
        23238195.5252303; ...
        21264390.0275339; ...
        21618588.1453816; ...
        20904050.8502602; ...
        23725281.3956461; ...
        24573685.5759472];
x_sv = [ 19053139.031138 , -6617959.60541346,  17102312.3210907; ...
        -4465692.00522436, -24635708.2535929,  8798897.30502178; ...
        -12135612.5430664, -23262031.6824513, -4239604.86332562; ...
         12382138.1241811, -22694316.604797 ,  5808363.86640313; ...
         8134457.44393974, -13576983.5497894,  21782090.2999914; ...
        -4740647.89490035, -15476222.1828813,  20881289.8519382; ...
        -18676899.6106716, -7432384.32812473,  17214937.0177376; ...
         2201167.23604536, -21921444.5956131, -14556724.9095099];
psr2 = [23223191.2941634; ...
        20670603.1445293; ...
        23238192.6971778; ...
        21264388.5397864; ...
        21618586.3513072; ...
        20904048.094015; ...
        23725277.84564; ...
        24573681.9139877];

% PART A - PERFECT CLOCK
fprintf("(a)\n");
x1 = [0;0;0];
e1 = 1;
x2 = [0;0;0];
e2 = 1;
while (e1 > 1e-6)
    % 2 sv solution
    u1 = x_sv(1:4,:) - x1';
    R1 = sqrt(sum(u1.^2,2));
    uv1 = u1./R1;

    H1 = -uv1;
    y1 = psr1(2:5) - R1;

    dx1 = pinv(H1)*y1;
    x1 = x1 + dx1;
    e1 = norm(dx1);
end
while (e2 > 1e-6)
    % 4 sv solution
    u2 = x_sv - x2';
    R2 = sqrt(sum(u2.^2,2));
    uv2 = u2./R2;

    H2 = -uv2;
    y2 = psr1 - R2;

    dx2 = pinv(H2)*y2;
    x2 = x2 + dx2;
    e2 = norm(dx2);
end
DOP1 = inv(H1'*H1);
DOP2 = inv(H2'*H2);
fprintf("lla1 = [%f, %f, %f]\n", ecef2lla(x1'));
fprintf("lla2 = [%f, %f, %f]\n", ecef2lla(x2'));
fprintf("PDOP1 = %f\n", sqrt(trace(DOP1)));
fprintf("PDOP2 = %f\n", sqrt(trace(DOP2)));

% PART B - IMPERFECT CLOCK
fprintf("(b)\n");
[x3,b3,~,DOP3,~] = gnssPos(x_sv(1:4,:), psr1(1:4), 1);
[x4,b4,~,DOP4,~] = gnssPos(x_sv, psr1, 1);
fprintf("lla3 = [%f, %f, %f]\n", ecef2lla(x3'));
fprintf("lla4 = [%f, %f, %f]\n", ecef2lla(x4'));
fprintf("PDOP3 = %f\n", sqrt(trace(DOP3(1:3,1:3))));
fprintf("PDOP4 = %f\n", sqrt(trace(DOP4(1:3,1:3))));

% PART C - SINGLE DIFFERENCE
fprintf("(c)\n");

% rcvr1 is reference
psr_ref = psr1;
x_ref = [x4;b4];
u_ref = x_sv - x_ref(1:3)';
R_ref = sqrt(sum(u_ref.^2,2));
uv_ref = u_ref ./ R_ref;

% rcvr2 is user
psr_user = psr2;
dR = psr_ref - psr_user;

H5 = [uv_ref(1:4,:), ones(4,1)];
y5 = dR(1:4);
dx5 = pinv(H5)*y5;
x5 = x_ref + dx5;
DOP5 = inv(H5'*H5);

H6 = [uv_ref, ones(8,1)];
y6 = dR;
dx6 = pinv(H6)*y6;
x6 = x_ref + dx6;
DOP6 = inv(H6'*H6);

fprintf("lla5 = [%f, %f, %f]\n", ecef2lla(x5(1:3)'));
fprintf("lla6 = [%f, %f, %f]\n", ecef2lla(x6(1:3)'));
fprintf("PDOP5 = %f\n", sqrt(trace(DOP5(1:3,1:3))));
% fprintf("GDOP5 = %f\n", sqrt(trace(DOP5)));
fprintf("PDOP6 = %f\n", sqrt(trace(DOP6(1:3,1:3))));
% fprintf("GDOP6 = %f\n", sqrt(trace(DOP6)));

% PART D - DOUBLE DIFFERENCE
fprintf("(d)\n");
H7 = uv_ref(2:4,:) - uv_ref(1,:);
y7 = dR(2:4) - dR(1);
dx7 = pinv(H7)*y7;
x7 = x_ref(1:3) + dx7;
DOP7 = inv(H7'*H7);

H8 = uv_ref(2:end,:) - uv_ref(1,:);
y8 = dR(2:end) - dR(1);
dx8 = pinv(H8)*y8;
x8 = x_ref(1:3) + dx8;
DOP8 = inv(H8'*H8);

fprintf("lla7 = [%f, %f, %f]\n", ecef2lla(x7(1:3)'));
fprintf("lla8 = [%f, %f, %f]\n", ecef2lla(x8(1:3)'));
fprintf("PDOP7 = %f\n", sqrt(trace(DOP7(1:3,1:3))));
fprintf("PDOP8 = %f\n", sqrt(trace(DOP8(1:3,1:3))));

tl = tiledlayout(2,2, Parent=tab(4), TileSpacing="tight");
title(tl, "\textbf{Geoplot of SV Positioning Tecniques}", Interpreter="latex", FontSize=18);
axi1 = geoaxes(Parent=tl);
lla1 = ecef2lla(x1');
lla2 = ecef2lla(x2');
lla3 = ecef2lla(x3');
lla4 = ecef2lla(x4');
lla5 = ecef2lla(x5(1:3)');
lla6 = ecef2lla(x6(1:3)');
lla7 = ecef2lla(x7');
lla8 = ecef2lla(x8');
hold(axi1, "on");
geoplot(lla1(1), lla1(2), 'o', LineWidth=3);
geoplot(lla2(1), lla2(2), 'o', LineWidth=3);
title("(A) Perfect Clocks");
geobasemap satellite;
nexttile;
geoplot(lla3(1), lla3(2), 'o', LineWidth=3);
hold on;
geoplot(lla4(1), lla4(2), 'o', LineWidth=3);
title("(B) Imperfect Clocks");
geobasemap satellite;
geolimits([32.586806, 32.587167], [-85.494444, -85.493333]);
nexttile;
geoplot(lla5(1), lla5(2), 'o', LineWidth=3);
hold on;
geoplot(lla6(1), lla6(2), 'o', LineWidth=3);
geoplot(lla4(1), lla4(2), 'x', LineWidth=3);
title("(C) Single Difference");
geobasemap satellite;
geolimits([32.586806, 32.587167], [-85.494444, -85.493333]);
nexttile;
geoplot(lla7(1), lla7(2), 'o', LineWidth=3);
hold on;
geoplot(lla8(1), lla8(2), 'o', LineWidth=3);
geoplot(lla4(1), lla4(2), 'x', LineWidth=3);
title("(D) Double Difference");
geobasemap satellite;
geolimits([32.586806, 32.587167], [-85.494444, -85.493333]);
l = legend({"4 SV", "8 SV", "Reference"}, Orientation="Horizontal");
l.Layout.Tile = "north";


%% PROBLEM 6
fprintf("\n<strong>PROBLEM 6</strong> \n");

% CH.2 PROB 1a-b PRN#4 PRN#7
PRN4 = caGen(5,9);
PRN7 = caGen(1,8);
[PRN4_corr, shifts] = correlation(PRN4, PRN4);

ax(5) = axes(Parent=tab(5));

t = tiledlayout(2,2, TileSpacing="compact");
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];
ax1 = axes(t);
stairs(ax1,1:16,PRN4(1:16), Color="#A2142F", LineWidth=3);
ax1.Box = 'off';
xlim(ax1,[1 16]);
ylim(ax1,[-1,2]);
ax2 = axes(t);
ax2.Layout.Tile = 2;
stairs(ax2,1007:1023,PRN4(1007:1023), color="#0072BD", LineWidth=3);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[1007 1023]);
ylim(ax2,[-1,2]);
linkaxes([ax1 ax2], 'y');
ylabel(ax1, "Binary Value");
title(bgAx, "\textbf{PRN4 First and Last 16 Chips}", Interpreter="latex");

bgAx2 = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx2.Layout.Tile = 3;
bgAx2.Layout.TileSpan = [1 2];
ax1 = axes(t);
ax1.Layout.Tile = 3;
stairs(ax1,1:16,PRN7(1:16), Color="#A2142F", LineWidth=3);
ax1.Box = 'off';
xlim(ax1,[1 16]);
ylim(ax1,[-1,2]);
ax2 = axes(t);
ax2.Layout.Tile = 4;
stairs(ax2,1007:1023,PRN7(1007:1023), color="#0072BD", LineWidth=3);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[1007 1023]);
ylim(ax2,[-1,2]);
linkaxes([ax1 ax2], 'y');
ylabel(ax1, "Binary Value");
title(bgAx2, "\textbf{PRN7 First and Last 16 Chips}", Interpreter="latex");

set(findall(gcf,'-property','FontSize'),'FontSize',16);

ax(6) = axes(Parent=tab(6));
plot(shifts, PRN4_corr ,LineWidth=2.5);
xlabel("Shifts");
ylabel("Correlation");
xlim([-1023, 1023]);
title("\textbf{PRN4 Corrleation of Chips 1:1023 and 1024:2046}", Interpreter="latex");
ax(6).FontSize = 16;

%% PROBLEM 7
fprintf("\n<strong>PROBLEM 7</strong> \n");
PRN4_ = PRN4;
PRN4_(PRN4_ == 0) = -1;
PRN7_ = PRN7;
PRN7_(PRN7_ == 0) = -1;

% PART A
ax(7) = axes(Parent=tab(7));
subplot(2,1,1);
h = histogram(PRN4_);
text([-1 0 1]-0.04, h.BinCounts+20, string(h.BinCounts));
title("\textbf{PRN4 Histogram}", Interpreter="latex");
subplot(2,1,2);
histogram(PRN7_);
title("\textbf{PRN7 Histogram}", Interpreter="latex");
text([-1 0 1]-0.04, h.BinCounts+20, string(h.BinCounts));
set(findall(gcf,'-property','FontSize'),'FontSize',16);

% PART B
tl = tiledlayout(2,1, TileSpacing="tight", Parent=tab(8));
xlabel(tl, "Chips/Samples", FontSize=16);
ylabel(tl, "Magnitude [dB]", FontSize=16);
ax(8) = axes(Parent=tl);
plot(20*log10(abs(fftshift(fft(PRN7_)))), LineWidth=2.5);
ax(8).XAxis.Visible = 'off';
grid on;
% periodogram(PRN4_);
title("\textbf{PSD of PRN4}", Interpreter="latex");
a = nexttile;
plot(20*log10(abs(fftshift(fft(PRN7_)))), LineWidth=2.5);
linkaxes([ax(8), a], 'x');
grid on;
xlim([0,1023]);
% periodogram(PRN7_);
title("\textbf{PSD of PRN7}", Interpreter="latex");
set(findall(gcf,'-property','FontSize'),'FontSize',16);

% PART C and D
tl = tiledlayout(3,1, TileSpacing="tight", Parent=tab(9));
xlabel(tl, "Shifts", FontSize=16);
ylabel(tl, "Correlation", FontSize=16);
ax(9) = axes(Parent=tl);
[PRN4_corr, shifts] = correlation(PRN4_, PRN4_);
plot(shifts, PRN4_corr ,LineWidth=2.5);
title("\textbf{PRN4 Auto Correlation}", Interpreter="latex", FontSize=16);
ax(9).XAxis.Visible = "off";
grid on;
a = nexttile;
[PRN7_corr, shifts] = correlation(PRN7_, PRN7_);
plot(shifts, PRN7_corr ,LineWidth=2.5);
title("\textbf{PRN7 Auto Correlation}", Interpreter="latex", FontSize=16);
a.XAxis.Visible = "off";
grid on;
b = nexttile;
[PRN4_7_corr, shifts] = correlation(PRN4_, PRN7_);
plot(shifts, PRN4_7_corr ,LineWidth=2.5);
linkaxes([ax(9), a, b], 'xy');
grid on;
xlim([-1023, 1023]);
ylim([-0.1, 1.1]);
title("\textbf{Cross Correlation}", Interpreter="latex", FontSize=16);
set(findall(gcf,'-property','FontSize'),'FontSize',16);


%% PROBLEM 8
fprintf("\n<strong>PROBLEM 8</strong> \n");

% 1023 chip c/a code = 1 ms
prn4 = caGen(5,9);
prn7 = caGen(1,8);
prn4(prn4==0) = -1;
prn7(prn7==0) = -1;

% L1 carrier sin wave
f_L1 = 1575.42e6;
t_L1 = 1/f_L1;

% carrier sampling and prn upsampling (lower frequency => higher time)
s_rate = 10;
s_per_c = 1575.42 / 1.023;
fs = f_L1*s_rate; % = f_chip
dt = 1/fs;
t = dt : dt : 0.001;

% carrier and prn upsampling
carrier = sin((2*pi*f_L1).*t);
prn4_ = repmat(repelem(prn4, int32(s_rate*s_per_c)),1);
prn7_ = repmat(repelem(prn7, int32(s_rate*s_per_c)),1);

% multiply
s4 = carrier.*prn4_;
s7 = carrier.*prn7_;

psd4 = 20*log10(abs(fftshift(fft(s4))));
psd7 = 20*log10(abs(fftshift(fft(s7))));

% plot
x_ax = linspace(f_L1-fs/2, f_L1+fs/2, length(carrier)) ./ 1e6 - 1575.42;

tl = tiledlayout(2,1, TileSpacing="tight", Parent=tab(10));
title(tl, "\textbf{PSD of Code Carrier Modulation}", Interpreter="latex", FontSize=18);
ylabel(tl, "Power [dB]", FontSize=16);
xlabel(tl, "Frequency [MHz]", FontSize=16)
ax00 = axes(Parent=tl);
plot(x_ax, psd4)
title("PRN4")
ax00.XAxis.Visible = "off";
ax00.FontSize=16;
grid on;
ax01 = nexttile;
plot(x_ax,psd7);
grid on;
title("PRN7")
ax01.FontSize=16;
linkaxes([ax00, ax01], 'xy');
xlim([f_L1/1e6-50, f_L1/1e6+50]);
ylim([0, 120])

% figure()
% periodogram(s4);

%%
% exportgraphics(tab(1), "./media/1b.png")
% exportgraphics(tab(2), "./media/1c.png")
% % exportgraphics(tab(3), "./media/4.png")
% exportgraphics(tab(4), "./media/5.png")
% exportgraphics(tab(5), "./media/6a.png")
% exportgraphics(tab(6), "./media/6b.png")
% exportgraphics(tab(7), "./media/7a.png")
% exportgraphics(tab(8), "./media/7b.png")
% exportgraphics(tab(9), "./media/7c.png")
% exportgraphics(tab(10), "./media/8.png")


