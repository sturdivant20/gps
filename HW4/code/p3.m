%% GPS HW3 - Problem 3 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 3</strong>\n");

% bevly encrypted signal
signal = generate_signal(2);
L = length(signal);

% controller design
BW = 1;                 % controller bandwidth (10-18 Hz)
K = 4;                  % loop gain
z = sind(45);           % damping ratio
w = BW*4*z / (1+z^2);   % natural frequency

% PI
Ki = w^2;
Kp = z*w;

% time setup
t_s = 1e-6;             % sampling time
t_int = 0.01;           % integration time
t = 0:t_s:t_int-t_s;    % simulated signal time
numSamp = t_int / t_s;  % samples per integrastion period
loop = 1:numSamp:L;

% initlize values
disc = zeros(size(loop));   % error
disc_int = 0;               % integrated error
eHat = zeros(size(loop));   % error estimate
theta = zeros(size(loop));  % phase angle
f = zeros(size(loop));      % frequency
I = zeros(size(loop));      % Inphase
Q = zeros(size(loop));      % Quadrature
f(1) = 100;

beg = 1;
fin = numSamp;

bit = zeros(1,88);
j = 0;

% figure;
for i = 2:length(loop)

    % data signal
    sig = signal(beg:fin);

    % Inphase and Quadrature
    I(i) = sig * sin(2*pi*f(i-1).*t' + theta(i-1));
    Q(i) = sig * cos(2*pi*f(i-1).*t' + theta(i-1));

    % error
    disc(i) = atan(Q(i) / I(i));
%     disc_int = disc_int + disc(i)*t_int;

    % propagate
%     f(i) = f(i-1) + Kp*disc(i) + Ki*disc_int;
    f(i) = f(i-1) + K * (Ki*t_int*disc(i) + Kp*(disc(i) - disc(i-1)));
    theta(i) = rem(theta(i-1) + 2*pi*f(i-1)*t_int, 2*pi);

    beg = beg + numSamp;
    fin = fin + numSamp;

%     plot(sig, 'r')
%     hold on;
%     plot(ssig, 'b')
%     legend("Data", "Sin")
%     hold off;
%     exportgraphics(gcf, 'PLL.gif', Append=true)


    % determine bit value
    if ~mod(i,100)
        j = j+1;
        if I(i) > 0
            bit(j) = 1;
        else
            bit(j) = 0;
        end
        
        % print character after 8 bits
        if ~mod(j,8) && j > 1
            bit_str = sprintf('%d', bit(j-7:j));
            letter = char(bin2dec(bit_str));
            fprintf("%c", letter);
        end
    end


end

figure
stairs(bit);
ylim([-0.2;1.2]);

% bit_str = sprintf('%d', bit);
% 
% fprintf("Message: ");
% for i = 1:8:length(bit)
%     letter = char(bin2dec(bit_str(i:i+7)));
%     fprintf("%c", letter);
% end
fprintf("\n");
