%% ECEF Tightly Coupled Extended Kalman Filter | Daniel Sturdivant
function [x,P] = ekfLoose(wb, fb, C, pos, vel, gps_pos, gps_vel, x, P, Q, R, dt, mode)
% ekfTight - Loosely coupled, error state, gnss-ins EKF
%
% Inputs:
%   wb      3x1     IMU angular velocity in BODY frame [rad/s]
%   fb      3x1     IMU specific force in BODY frame [m/s^2]
%   C       3x3     INS Rotation Matrix
%   pos     3x1     INS user position in ECEF [m]
%   vel     3x1     INS user velocity in ECEF [m]
%   psr     Mx1     GPS psuedorange measurements [m]
%   dopp    Mx1     GPS psuedorange-rate measurements [m/s]
%   x_sv    Mx3     GPS satellite positions in ECEF [m]
%   x_sv    Mx3     GPS satellite velocities in ECEF [m/s]
%   x       17x1    Current KF Error Sate (x+(k-1))
%   P       17x17   Current KF Covariance (P+(k-1))
%   Q       17x17   Process Noise
%   R       MxM     Measurement Noise
%   dt      1x1     Time step [s]
%   mode    str     'time' or 'meas' or 'all'
%
% Output:
%   x       17x1    Updated state (x+(k))
%   P       17x1    Update covariance (P+(k))
%

%% TIME UPDATE
if strcmp(mode, 'all') || strcmp(mode, 'time')

    % rotation rate of earth in ecef
    W_ie_e = utils.vec2skew([0; 0; 7.292115e-5]);

    % geocentric radius (Groves 2.137)
    lla = ecef2lla(pos');
    L = deg2rad(lla(1));
    e2 = 0.00669438000426081;
    [~, Re] = utils.radii(L);
    r_eS_e = Re * sqrt(cos(L)^2 + (1 - e2)^2 * sin(L)^2);
    
    % gravity
    g = utils.gravityECEF(pos);
    
    % state transition matrix (Groves 14.49)
    A = zeros(15,15);
    A(1:3,1:3) = -W_ie_e;
    A(1:3,13:15) = C;
    A(4:6,1:3) = -utils.vec2skew(C*fb);
    A(4:6,4:6) = -2*W_ie_e;
    A(4:6,7:9) = -(2 * g / r_eS_e) * (pos' / sqrt(pos' * pos));
    A(4:6,10:12) = C;
    A(7:9,4:6) = eye(3);

    % control matrix (dr. aly el-osery -> aided INS)
%     B = zeros(15,15);
%     B(1:3,1:3) = -C;
%     B(4:6,4:6) = -C;
%     B(10:15,10:15) = eye(6);

    % discretization
    A = eye(15) + A*dt;
%     Q = (A*B*Q*B'*A' + B*Q*B') * dt/2;

    % Propagation
    x = zeros(15,1);
    P = A*P*A' + Q;

end


%% MEASUREMENT UPDATE
if strcmp(mode, 'all') || strcmp(mode, 'meas')

    % generate predicted measurements
    y = [gps_pos - pos; ...
         gps_vel - vel];

    Z = zeros(3);
    I = eye(3);

    % nonlinear measurement matrix (Groves 14.126)
    C = [Z,  Z, -I, Z, Z; ...
         Z, -I,  Z, Z, Z];

    % Correction
    L = P*C'*(C*P*C' + R)^-1;
    P = (eye(15) - L*C)*P;
    x = L*y;

end

end