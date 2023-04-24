%% GPS HW3 - Problem 5 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 5</strong>\n");


% METHOD
% method = 'serial';
method = 'parallel';


% data information
f_s = 20e6;  % 20 MHz [Hz]
T_s = 0.001; % 1 ms of data [s]
f_if = 5000445.88565834;
nSamples = f_s*T_s;
n_ms = 10;

% read binary file
fName = "gpsBase_IFEN_IF_new2.bin";
fid = fopen(fName);
fseek(fid, 0, 'bof');

% load in PRN
load("ca_codes.mat");

%%

% upsampling indices
idx_ = 1/nSamples : 1023/nSamples : 1023-(1/nSamples);
idx = ceil(idx_);
idx = repmat(idx, 1, n_ms);
t_samp = (0 : 1/nSamples*T_s : n_ms*T_s-(1/nSamples*T_s))';
t_prn = 0 : 1/1023*T_s : T_s - (1/1023*T_s);

% % % upsampling
% prn = ca_code(7,:)';
% % prn_up = zeros(20000,1);
% % ttt = 1;
% % for tt = 1:20000
% %     if (t_samp(tt) >= t_prn(ttt)) && (t_samp(tt) < t_prn(ttt+1))
% %         prn_up(tt) = prn(ttt);
% %     else
% %         ttt = ttt + 1;
% %         if ttt == 1024
% %             ttt = 1;
% %         end
% %         prn_up(tt) = prn(ttt);
% %     end
% % end
% % prn_up = circshift(prn_up,20);
% 
% f = figure;
% tbs = uitabgroup(Parent=f);
% tab(1) = uitab(Parent=tbs, Title="PRN");
% tab(2) = uitab(Parent=tbs, Title="Xcorr");
% axes(Parent=tab(1));
% hold on
% stairs(t_prn, prn)
% % stairs(t_samp, prn_up, '*')
% stairs(t_samp, prn(idx), '*')
% ylim([-1.5,1.5])
% legend("PRN", "Upsampled")
% axes(Parent=tab(2));
% % plot(xcorr(prn_up))
% plot(xcorr(prn(idx)))


%%
prn_found = logical(zeros(1,32));
known_prn = [7,19,30];

f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
b = 1;


if strcmp(method, 'serial')

% initial doppler offset
f_dopp = -10000:500:10000;

% auto correlation
R1 = zeros(1023,length(f_dopp),32);

% NUMBER OF DATA SAMPLES TO USE
for a = 1:10

    % n_ms ms of data
    [sig1, ~] = fread(fid, n_ms*nSamples, 'int8');

    % PRNS IN SEARCH
    for k = known_prn
        if prn_found(k) == 0
    
            % upsampling
            prn = [ca_code(1,1)'; ca_code(1,:)'];
            prn_up = zeros(20000,1);
            ttt = 1;
            for tt = 1:20000
                if t_samp(tt) <= t_prn(ttt)
                    prn_up(tt) = prn(ttt);
                else
                    ttt = ttt + 1;
                    prn_up(tt) = prn(ttt);
                end
            end
        
            % PRN SHIFT
            for i = 1:1023
            
                tau = i-1;
                prn_test = circshift(prn_up, tau*20);
            
                % DOPPLER BINS
                for j = 1:length(f_dopp)
            
                    % what bevly said
                    f = f_if + f_dopp(j);
                    I = sig1 .* prn_test .* sin(2*pi*f*t);
                    Q = sig1 .* prn_test .* cos(2*pi*f*t);
                    R1(i,j,k) = sum(I)^2 + sum(Q)^2;
            
                end
            
            end
        
            % detection
            tmp1 = sort(max(R1(:,:,k), [], 1), 'descend');
            ratio1 = tmp1(1) / tmp1(2);
        
            % ratio test
            if ratio1 > 5
                fprintf("SV #%d acquired!\n", k);
                prn_found(k) = 1;

                ttl(b) = sprintf("SV%d", k);
                tab(b) = uitab(Parent=tbs, Title=ttl(b));
                ax1 = axes(Parent=tab(b));
                surf(ax1, f_dopp, 0:1022, R1(:,:,k));
                title(sprintf("SV %d Acquisition", k));
                xlabel('Doppler Shift [kHz]');
                ylabel('PRN Shift [chips]');
                zlabel('Correlation')

                b = b+1;
            end

        end
    end

end

elseif strcmp(method, 'parallel')

% initial doppler offset
f_dopp = -10000:100:10000;

% auto correlation
R1 = zeros(n_ms*nSamples,length(f_dopp),32);

% NUMBER OF DATA SAMPLES TO USE
for a = 1:1

    % n_ms ms of data
    [sig1, ~] = fread(fid, n_ms*nSamples, 'int8');

    % PRNS IN SEARCH
    for k = 7:7
        if prn_found(k) == 0

            % upsampling
            prn = ca_code(k,:)';
            prn_up = prn(idx);

            % DOPPLER BINS
            for j = 1:length(f_dopp)
        
                % what bevly said
                f = f_if + f_dopp(j);
                I = sig1 .* sin(2*pi*f*t_samp);
                Q = sig1 .* cos(2*pi*f*t_samp);
                R1(:,j,k) = abs(ifft(fft(I + Q.*1j) .* conj(fft(prn_up)))).^2;
        
            end

            % detection
            tmp1 = sort(max(R1(:,:,k), [], 1), 'descend');
            ratio1 = tmp1(1) / tmp1(2);
        
            % ratio test
%             if ratio1 > 5
                fprintf("SV #%d acquired!\n", k);
                prn_found(k) = 1;

                ttl(b) = sprintf("SV%d", k);
                tab(b) = uitab(Parent=tbs, Title=ttl(b));
                ax1 = axes(Parent=tab(b));
                mesh(ax1, f_dopp/1e3, idx_, R1(1:nSamples,:,k));
                title(sprintf("SV %d Acquisition", k));
                xlabel('Doppler Shift [kHz]');
                ylabel('PRN Shift [chips]');
                zlabel('Correlation')

                b = b+1;
%             end

        end
    end

end

end

%% 
for i = 1:b-1
    exportgraphics(tab(i), sprintf('./media/p7_%s_2.png', ttl(i)));
end

