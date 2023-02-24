function [clkCorr, pos, vel] = svPVT(ephemeris, transmitTime, transitTime)
%
% Inputs:
%   ephemeris   Ephemeris of all satellites in use
%   transitTime Range/Speed of light of all satellites in use
%
% Outputs:
%   clkCorr     clock correction
%   pos         Satellite positions in ECEF
%   vel         Satellite velocities in ECEF
%

% load('RCVR_S1_data.mat')
% 
% remove NaN
% ephemeris = cell2mat(struct2cell(RCVR_S1.ephem))';
% cond = ~isnan(ephemeris(:,2)); % condition for removing NaN in array
% ephemeris = ephemeris(cond, :);

% constants
GM = 3.986005e14;             % Earth's universal gravitation parameter
OmegaDot_E = 7.2921151467e-5;   % Semi-Major axis rotation rate
c = 299792458;                  % speed of light
Pi = 3.1415926535898;           % GPS Pi
F = -4.442807633e-10;           % ????????? [sec/m^(1/2)]

% ephemeris variables
prn = ephemeris(:,1);
a_f0 = ephemeris(:,2);
a_f1 = ephemeris(:,3);
a_f2 = ephemeris(:,4);
IODE_sf2 = ephemeris(:,5);
C_rs = ephemeris(:,6);
deltan = ephemeris(:,7);
M_0 = ephemeris(:,8);
C_uc = ephemeris(:,9);
e = ephemeris(:,10);
C_us = ephemeris(:,11);
A = ephemeris(:,12).^2;
t_oe = ephemeris(:,13);
C_ic = ephemeris(:,14);
Omega_0 = ephemeris(:,15);
C_is = ephemeris(:,16);
i_0 = ephemeris(:,17);
C_rc = ephemeris(:,18);
omega = ephemeris(:,19);
OmegaDot = ephemeris(:,20);
IDOT = ephemeris(:,21);
L2Code = ephemeris(:,22);
gpsWeek = ephemeris(:,23);
L2P = ephemeris(:,24);
accuracy = ephemeris(:,25);
health = ephemeris(:,26);
T_GD = ephemeris(:,27);
IODC = ephemeris(:,28);
t_oc = ephemeris(:,29);

% CLOCK TOW ADJUSTMENT
dt = time_check(transmitTime, t_oc);
clkCorr = a_f0 + a_f1.*dt + a_f2.*dt.^2 - T_GD;
tk = time_check(transmitTime - clkCorr, t_oe);

% MOTION
n_0 = sqrt(GM ./ A.^3); % computed mean motion [rad/s]
n = n_0 + deltan;       % corrected mean motion

% ANOMOLY
M = rem((M_0 + n.*tk) + 2*Pi, 2*Pi);% Mean anomoly
error = Inf;
E_old = M;                          % Eccentric anomaly
while all(error > 1e-12)
    E = E_old + ((M - E_old +  e.*sin(E_old)) ./ (1 - e.*cos(E_old)));
%     E = M + e.*sin(E_old);
    error = abs(rem(E - E_old, 2*Pi));
    E_old = E;
end
E = rem(E + 2*Pi, 2*Pi);
nu = 2.*atan( sqrt((1+e)./(1-e)) .* tan(E./2) );    % true anamoly
% nu = atan2(sqrt(1 - e.^2) .* sin(E), cos(E) - e);

% SV POSITION
Phi = rem(nu + omega, 2*Pi);                % argument of latitiude
du = C_us.*sin(2.*Phi) + C_uc.*cos(2.*Phi); % argument of latitude correction
dr = C_rs.*sin(2.*Phi) + C_rc.*cos(2.*Phi); % radius correction
di = C_is.*sin(2.*Phi) + C_ic.*cos(2.*Phi); % inclination correction
u = Phi + du;                               % corrected argument of latitude
r = A.*(1 - e.*cos(E)) + dr;                % corrected radius
i = i_0 + di + IDOT.*tk;                    % corrected inclination
x_ = r.*cos(u);                             % x in orbital frame
y_ = r.*sin(u);                             % y in orbital frame
Omega_ = Omega_0 + (OmegaDot - OmegaDot_E).*tk - OmegaDot_E.*(t_oe + transitTime); % corrected longitude of ascending node
Omega_ = rem(Omega_ + 2*Pi, 2*Pi);

% ecef position
x = x_.*cos(Omega_) - y_.*cos(i).*sin(Omega_);
y = x_.*sin(Omega_) + y_.*cos(i).*cos(Omega_);
z = y_.*sin(i);
pos = [x, y, z];

% SV VELOCITY
EDot = n ./ (1-e.*cos(E));  % Eccentric anomoly rate
nuDot = (EDot.*sqrt(1 - e.^2)) ./ (1 - e.*cos(E));  % true anomly rate
di_dt = IDOT + 2.*nuDot.*(C_is.*cos(2.*Phi) - C_ic.*sin(2.*Phi));   % corrected inclination angle rate
uDot = nuDot + 2.*nuDot.*(C_rs.*cos(2.*Phi) - C_uc.*sin(2.*Phi));   % corrected argument of latitude rate
rDot = e.*A.*EDot.*sin(E) + 2.*nuDot.*(C_rs.*cos(2.*Phi) - C_rc.*sin(2.*Phi)); % corrected radius rate

OmegaDot_ = OmegaDot - OmegaDot_E;  % longitude of ascending node rate
xDot_ = rDot.*cos(u) - r.*uDot.*sin(u); % orbital x velocity
yDot_ = rDot.*sin(u) + r.*uDot.*cos(u); % orbital y velocity

% ecef velocity
xDot = -x_.*OmegaDot_.*sin(Omega_) + xDot_.*cos(Omega_) - yDot_.*sin(Omega_).*cos(i) - y_.*(OmegaDot_.*cos(Omega_).*cos(i) - di_dt.*sin(Omega_).*sin(i));
yDot = x_.*OmegaDot_.*cos(Omega_) + xDot_.*sin(Omega_) + yDot_.*cos(Omega_).*cos(i) - y_.*(OmegaDot_.*sin(Omega_).*cos(i) + di_dt.*cos(Omega_).*sin(i));
zDot = yDot_.*sin(i) + y_.*di_dt.*cos(i);
vel = [xDot, yDot, zDot];

% CLOCK CORRECTIONS
dtr = F .* e .* sqrt(A) .* sin(E);
clkCorr = a_f0 + a_f1.*dt + a_f2.*dt.^2 + dtr - T_GD;


% correct for ToW crossover
function tk = time_check(t, t_oe)
    % time of ephemeris reference epoch
    tk = t - t_oe;

    % correct time
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
end

end
