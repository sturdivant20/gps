%% ECEF Tightly Coupled Extended Kalman Filter | Daniel Sturdivant
function [x,P] = ekfTight(wb, fb, C, pos, vel, psr, dopp, x_sv, v_sv, x, P, Q, R, dt, mode)
% ekfTight - Tightly coupled, error state, gnss-ins EKF
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


% % rotation rate of earth in ecef
w_ie_e = 7.292115e-5;
W_ie_e = utils.vec2skew([0;0;w_ie_e]);

%% TIME UPDATE
if strcmp(mode, 'all') || strcmp(mode, 'time')

    % geocentric radius (Groves 2.137)
    lla = ecef2lla(pos');
    L = deg2rad(lla(1));
    e2 = 0.00669438000426081;
    [~, Re] = utils.radii(L);
    r_eS_e = Re * sqrt(cos(L)^2 + (1 - e2)^2 * sin(L)^2);
    
    % gravity
    g = utils.gravityECEF(pos);
    
    % state transition matrix (Groves 14.49)
    A = zeros(17,17);
    A(1:3,1:3) = -W_ie_e;
    A(1:3,13:15) = C;
    A(4:6,1:3) = -utils.vec2skew(C*fb);
    A(4:6,4:6) = -2*W_ie_e;
    A(4:6,7:9) = -(2 * g / r_eS_e) * (pos' / sqrt(pos' * pos));
    A(4:6,10:12) = C;
    A(7:9,4:6) = eye(3);
    A(16,17) = 1;

    % control matrix (dr. aly el-osery -> aided INS)
%     B = zeros(17,17);
%     B(1:3,1:3) = -C;
%     B(4:6,4:6) = -C;
%     B(10:17,10:17) = eye(8);

    % discretization
    A = eye(17) + A*dt;
%     Q = (A*B*Q*B'*A' + B*Q*B') * dt/2;

    % Propagation
    x = [zeros(15,1); x(16)+x(17)*dt; x(17)];
    P = A*P*A' + Q;

end


%% MEASUREMENT UPDATE
if strcmp(mode, 'all') || strcmp(mode, 'meas')

    % generate predicted measurements
    c = 299792458;          % speed of light [m/s]
    len = length(psr);
    u = zeros(len,3);
    yHat = zeros(2*len,1);

    for i = 1:len
%         % approximate range
%         dr = x_sv(i,:)' - pos;
%         r_approx = sqrt(dr' * dr);
% 
%         % signal transit rotation (Groves 8.36)
%         C_e_I = [                        1, w_ie_e(3) * r_approx / c, 0;...
%                  -w_ie_e(3) * r_approx / c,                        1, 0;...
%                                          0,                        0, 1];
% 
%         % predict psuedorange (Groves 9.165)
%         dr = C_e_I * x_sv(i,:)' - pos;
%         r = sqrt(dr' * dr);
% 
%         % unit vector
%         u(i,:) = dr' ./ r;
% 
%         % predict doppler (Groves 9.165)
%         v = u(i,:) * (C_e_I * (v_sv(i,:)' + W_ie_e*x_sv(i,:)') - (vel + W_ie_e*pos));

        dr = x_sv(i,:)' - pos;
        r = sqrt(dr' * dr);
        u(i,:) = dr' ./ r;
        v = u(i,:) * (v_sv(i,:)' - vel);

        % estimated measurement
        yHat(i) = r + x(16);
        yHat(i+len) = v + x(17);
    end

    Z = zeros(len,3);
    z = zeros(len,1);
    o = ones(len,1);

    % nonlinear measurement matrix (Groves 14.126)
    C = [Z, Z, u, Z, Z, o, z; ...
         Z, u, Z, Z, Z, z, o];
    y = [psr; dopp];

    % Correction
    L = P*C'*(C*P*C' + R)^-1;
    P = (eye(17) - L*C)*P;
    x = x + L*(y - yHat);

end

end