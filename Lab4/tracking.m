%% INIT
clc; close all; clear;
fprintf("<strong>INIT</strong>\n");

% data information
f_s = 20e6;  % 20 MHz [Hz]
T_s = 0.001; % 1 ms of data [s]
t_int = T_s;
f_if = 5000445.88565834;
nSamples = f_s*T_s;
n_ms = 1;

% read binary file
fName = 'gpsBase_IFEN_IF_new2.bin';
fid = fopen(sprintf('%s',fName));
fseek(fid, 0, 'bof');

% load in PRN 7
load("ca_codes.mat");

fprintf("FileName: '%s' \n", fName);
fprintf("Sample Rate: %de+6 Hz \n", f_s*1e-6);
fprintf("Integration Time: %de-3 s \n", t_int*1e3);
fprintf("\n");


%% ACQUISITION
fprintf("<strong>ACQUISITION</strong>\n");

svInUse = [1,7,14,17,21,30];

% n_ms ms of data
n_ms = 1;
[sig, ~] = fread(fid, n_ms*20000, 'int8');

% time for genrated signal
t = 0 : 1/nSamples*T_s : n_ms*T_s-(1/nSamples*T_s);

f_dopp = -5000:10:5000;
R = zeros(length(f_dopp), n_ms*nSamples, length(svInUse));

% initialize (for speed :])
dopp_offset = zeros(1,length(svInUse));
code_offset = zeros(1,length(svInUse));

k = 1;
for i = svInUse
    % upsampled prn
    [prn_up, ~, ~] = upsampleCode(ca_code(i,:), 20e6, 1.023e6, 0);
    prn_up = repmat(prn_up,1,n_ms);

    % DOPPLER BINS
    for j = 1:length(f_dopp)
    
        f_phase = f_if + f_dopp(j);
        I = sig' .* sin(2*pi*f_phase*t);
        Q = sig' .* cos(2*pi*f_phase*t);
        R(j,:,k) = abs(ifft( fft(I + Q.*1j) .* conj(fft(prn_up)) )).^2;
    
    end

    % initial frequency guess
    [dopp_offset(k),code_offset(k),~] = find(R(:,:,k) == max(R(:,1:20000,k),[],'all'),1);
    fprintf("SV%d\n",i);
    fprintf("Initial Doppler: %d Hz\n", f_dopp(dopp_offset(k)));
    fprintf("Initial Code Offset: %d chips\n", round(code_offset(k)/20));
    fprintf("---------------------------------------------------------\n");
    k = k+1;

end

% f = figure;
% tbs = uitabgroup(Parent=f);
% for i = 1:length(svInUse)
%     tab(i) = uitab(Parent=tbs, Title=sprintf("SV%d", svInUse(i)));
%     axes(Parent=tab(i));
%     mesh(linspace(1,1023,20000), f_dopp/1e3, R(:,1:20000,i));
%     title(sprintf("SV %d Acquisition", svInUse(i)));
%     ylabel('Doppler Shift [kHz]');
%     xlabel('PRN Shift [chips]');
%     zlabel('Correlation')
% end

clearvars -except fid sig f_dopp dopp_offset code_offset svInUse f_if f_s ca_code
fprintf("\n");


%% TRACKING
fprintf("<strong>TRACKING</strong>\n");

fName = 'gpsBase_IFEN_IF_new2.bin';
fid = fopen(sprintf('%s',fName));
fseek(fid, 0, 'bof');

% PLL
BW = 25;                   % PLL bandwidth (10-18 Hz)
K1 = 4;                    % loop gain
z = 1.0*sind(45);          % damping ratio
w = BW*4*z / (1+z^2);      % natural frequency
Kp_pll = z*w;
Ki_pll = w^2;

% DLL
BW = 10;                    % DLL bandwidth (1-5 Hz)
K2 = 4;                    % loop gain
z = 2.0*sind(45);          % damping ratio
w = BW*4*z / (1+z^2);      % natural frequency
Kp_dll = z*w;
Ki_dll = w^2;

% initialize (for speed :])
len = 1*1000;      % 60 seconds
IE = zeros(length(svInUse),len);
IP = zeros(length(svInUse),len);
IL = zeros(length(svInUse),len);
QE = zeros(length(svInUse),len);
QP = zeros(length(svInUse),len);
QL = zeros(length(svInUse),len);
f_chip = zeros(length(svInUse),len);
f_phase = zeros(length(svInUse),len);
theta = zeros(length(svInUse),len);
codePhase = zeros(length(svInUse),len);
dataSize = zeros(length(svInUse),len);
bit = nan(length(svInUse), len/20);

f_chip(:,1) = 1.023e6 + 1.023e6/1575.42e6*f_dopp(dopp_offset)';
f_phase(:,1) = f_if + f_dopp(dopp_offset)';

k = zeros(length(svInUse),1);       % index where first preamble occurs [ms]
ii = zeros(length(svInUse),1);      % bit index [20 ms]
bit_new = zeros(length(svInUse),1);
frame = logical(zeros(length(svInUse),1));
disc_pll = zeros(length(svInUse),1);
disc_dll = zeros(length(svInUse),1);
preamble = [1, -1, -1, -1, 1, -1, 1, 1];

for j = 2

    for i = 1:len

        % upsample prn and read in data
        [prn_up, codePhase(j,i+1), dataSize(j,i)] = upsampleCode(ca_code(svInUse(j),:), f_s, f_chip(j,i), codePhase(j,i));
        [sig, ~] = fread(fid, dataSize(j,i), 'int8');
    
        t_int = dataSize(j,i) / 20e6;
        t = (0 : 1/dataSize(j,i) : 1 - 1/dataSize(j,i)) * t_int;
    
        % inphase and quadrature
        IE(j,i) = sum(circshift(prn_up, code_offset+10) .* sig' .* sin(2*pi*f_phase(j,i)*t + theta(j,i)));
        IP(j,i) = sum(circshift(prn_up, code_offset)    .* sig' .* sin(2*pi*f_phase(j,i)*t + theta(j,i)));
        IL(j,i) = sum(circshift(prn_up, code_offset-10) .* sig' .* sin(2*pi*f_phase(j,i)*t + theta(j,i)));
        QE(j,i) = sum(circshift(prn_up, code_offset+10) .* sig' .* cos(2*pi*f_phase(j,i)*t + theta(j,i)));
        QP(j,i) = sum(circshift(prn_up, code_offset)    .* sig' .* cos(2*pi*f_phase(j,i)*t + theta(j,i)));
        QL(j,i) = sum(circshift(prn_up, code_offset-10) .* sig' .* cos(2*pi*f_phase(j,i)*t + theta(j,i)));
    
        % dll
        disc_dll_old = disc_dll(j);
        L = sqrt(IL(j,i)^2 + QL(j,i)^2);
        E = sqrt(IE(j,i)^2 + QE(j,i)^2);
        disc_dll(j) = 0.5 * (E-L) / (E+L);
        f_chip(j,i+1) = f_chip(j,i) - K2*(Ki_dll*t_int*disc_dll(j) + Kp_dll*(disc_dll(j) - disc_dll_old));
    
        % pll
        disc_pll_old = disc_pll(j);
        disc_pll(j) = atan(QP(j,i) / IP(j,i));
        theta(j,i+1) = rem(theta(j,i) + 2*pi*f_phase(j,i)*t_int, 2*pi);
        f_phase(j,i+1) = f_phase(j,i) + K1*(Ki_pll*t_int*disc_pll(j) + Kp_pll*(disc_pll(j) - disc_pll_old));
    
%         % record bit
%         if k(j) == 0
%             bit_old = bit_new(j);
%             bit_new(j) = sign(IP(j,i));
%         
%             % determine first bit flip
%             if i > 500
%                 if bit_old ~= bit_new(j)
%                     fprintf("SV%d First real bit flip: %d ms\n", svInUse(j), i);
%                     k(j) = i;
%                     bit(j,1) = bit_new(j);
%                     ii(j) = 2;
%                 end
%             end
%         else
%             if ~mod(i,20)
%                 bit(j,ii(j)) = sign(IP(j,i-1));
%     
%                 % find preamble back in time
%                 if ii(j) > 300 && ~frame(j)
%                     if isequal(bit(j,ii(j)-7:ii(j)), preamble) || isequal(bit(j,ii(j)-7:ii(j)), -preamble)
%                         if isequal(bit(j,(ii(j)-7:ii(j))-300), preamble) || isequal(bit(j,(ii(j)-7:ii(j))-300), -preamble)
%                             fprintf("SV%d Preamble starts at bit %d\n", svInUse(j), ii(j)-307);
%                             frame(j) = true;
%                             k(j) = ii(j)-307;
%                         end
%                     end
%                 end
%     
%                 ii(j) = ii(j) + 1;
%             end
%         end

    end

end

% close file
fclose(fid);
fprintf("\n");


% %% data
% subframes?
% b1.frame = bit(258:258+299);
% b1.frame = bit(j,k(j):k(j)+299);
% b1.frame(b1.frame < 0) = 0;
% b1.preamble = b1.frame(1:8);
% b1.Fnum = b1.frame(50:52);
% b1.ToW = b1.frame(31:47);
% b1.Fnum_n = bin2dec(sprintf('%d', b1.Fnum));
% b1.ToW_n = 6*bin2dec(sprintf('%d', b1.ToW));
% 
% b2.frame = bit(558:558+299);
% b2.frame = bit(j,(k(j):k(j)+299)+300);
% b2.frame(b2.frame < 0) = 0;
% b2.preamble = b2.frame(1:8);
% b2.Fnum = b2.frame(50:52);
% b2.ToW = b2.frame(31:47);
% b2.Fnum_n = bin2dec(sprintf('%d', b2.Fnum));
% b2.ToW_n = 6*bin2dec(sprintf('%d', b2.ToW));
% 
% b3.frame = bit(858:858+299);
% b3.frame = bit(j,(k(j):k(j)+299)+600);
% b3.frame(b3.frame < 0) = 0;
% b3.preamble = b3.frame(1:8);
% b3.Fnum = b3.frame(50:52);
% b3.ToW = b3.frame(31:47);
% b3.Fnum_n = bin2dec(sprintf('%d', b3.Fnum));
% b3.ToW_n = 6*bin2dec(sprintf('%d', b3.ToW));
% 
% b4.frame = bit(1158:1158+299);
% b4.frame = bit(j,(k(j):k(j)+299)+900);
% b4.frame(b4.frame < 0) = 0;
% b4.preamble = b4.frame(1:8);
% b4.Fnum = b4.frame(50:52);
% b4.ToW = b4.frame(31:47);
% b4.Fnum_n = bin2dec(sprintf('%d', b4.Fnum));
% b4.ToW_n = 6*bin2dec(sprintf('%d', b4.ToW));
% 
% b5.frame = bit(1458:1458+299);
% b5.frame = bit(j,(k(j):k(j)+299)+1200);
% b5.frame(b5.frame < 0) = 0;
% b5.preamble = b5.frame(1:8);
% b5.Fnum = b5.frame(50:52);
% b5.ToW = b5.frame(31:47);
% b5.Fnum_n = bin2dec(sprintf('%d', b5.Fnum));
% b5.ToW_n = 6*bin2dec(sprintf('%d', b5.ToW));
% 
% b6.frame = bit(1758:1758+299);
% b6.frame = bit(j,(k(j):k(j)+299)+1500);
% b6.frame(b6.frame < 0) = 0;
% b6.preamble = b6.frame(1:8);
% b6.Fnum = b6.frame(50:52);
% b6.ToW = b6.frame(31:47);
% b6.Fnum_n = bin2dec(sprintf('%d', b6.Fnum));
% b6.ToW_n = 6*bin2dec(sprintf('%d', b6.ToW));
% 
% b7.frame = bit(2058:2058+299);
% b7.frame = bit(j,(k(j):k(j)+299)+1800);
% b7.frame(b7.frame < 0) = 0;
% b7.preamble = b7.frame(1:8);
% b7.Fnum = b7.frame(50:52);
% b7.ToW = b7.frame(31:47);
% b7.Fnum_n = bin2dec(sprintf('%d', b7.Fnum));
% b7.ToW_n = 6*bin2dec(sprintf('%d', b7.ToW));
% 
% b8.frame = bit(2358:2358+299);
% b8.frame = bit(j,(k(j):k(j)+299)+2100);
% b8.frame(b8.frame < 0) = 0;
% b8.preamble = b8.frame(1:8);
% b8.Fnum = b8.frame(50:52);
% b8.ToW = b8.frame(31:47);
% b8.Fnum_n = bin2dec(sprintf('%d', b8.Fnum));
% b8.ToW_n = 6*bin2dec(sprintf('%d', b8.ToW));
% 
% b9.frame = bit(2658:2658+299);
% b9.frame = bit(j,(k(j):k(j)+299)+2400);
% b9.frame(b9.frame < 0) = 0;
% b9.preamble = b9.frame(1:8);
% b9.Fnum = b9.frame(50:52);
% b9.ToW = b9.frame(31:47);
% b9.Fnum_n = bin2dec(sprintf('%d', b9.Fnum));
% b9.ToW_n = 6*bin2dec(sprintf('%d', b9.ToW));


% %%
figure;
hold on;
plot(IP(j,:), '.', LineWidth=1.5);
plot(QP(j,:), '.', LineWidth=1.5);
title("Tracking Prompt Correlators");
xlabel("Integration Periods [1 ms]");
ylabel("Correlation");
legend("IP", "QP");

figure;
hold on;
plot(IE(j,:), '.', LineWidth=1.5);
plot(IL(j,:), '.', LineWidth=1.5);
plot(QE(j,:), '.', LineWidth=1.5);
plot(QL(j,:), '.', LineWidth=1.5);
title("Tracking Early / Late Correlators");
xlabel("Integration Periods [1 ms]");
ylabel("Correlation");
legend("IE", "IL", "QP", "QL");
% 
% figure;
% stairs(bit);
% title("Data Bits");
% ylim([-1.1,1.1]);
