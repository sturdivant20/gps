%% ECEF IMU Mechanization | Daniel Sturdivant
function [P,V,C] = imuMechECEF(w_ib_b, f_ib_b, P_old, V_old, C_old, dt)
% IMUMECH - IMU body-to-ecef mechanization
%
% Inputs:
%   w_ib_b  3x1     IMU angular velocty X,Y,Z [rad/s]
%   f_ib_b  3x1     IMU specific force / acceleration X,Y,Z [m/s^2]
%   P_old   3x1     Previous ECEF Position X,Y,Z [m]
%   V_old   3x1     Previous ECEF Velocity X,Y,Z [m/s]
%   C_old   3x3     Previous body-to-ecef rotation matrix
%   dt      1x1     time step [s]
%
% Outputs:
%   P       3x1     updated position [m]
%   V       3x1     updated velcoity [m/s]
%   A       3x1     updated attitude (euler angles) [rad]
%   c       3x3     updated rotation matrix
%   q       4x1     updated quaternion
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integr-ated Navigation 
%       Systems (2013)  -  Paul D. Groves
%

%% ATTITUDE UPDATE

% Groves 2.145
w_ie = 7.292115e-5;
W_ie_e = utils.vec2skew([0,0,w_ie]);
a_ie = w_ie * dt;
A_ie_e = utils.vec2skew([0,0,a_ie]);
C_earth = [ cos(a_ie), sin(a_ie), 0; ...
           -sin(a_ie), cos(a_ie), 0;
                    0,         0, 1];

% Groves 5.69/5.73 (Rodrigues' Formula)
Alpha = utils.vec2skew(w_ib_b * dt);
C_bar = expm(Alpha);

% Update attitude (Groves 5.75)
C = C_earth * C_old * C_bar;


%% VELOCITY UPDATE

% specific force transformation (Groves 5.83/5.84/5.85)
C_avg = C_old * expm(Alpha./2) - 0.5*A_ie_e*C_old;
f_ib_e = C_avg * f_ib_b;

% Velocity Update (Groves 5.36)
g = utils.gravityECEF(P_old);
% g = zeros(3,1);
% [g(1), g(2), g(3)] = gravitysphericalharmonic(P_old');
V = V_old + (f_ib_e + g - 2 * W_ie_e * V_old) * dt;


%% POSITION UPDATE

% Groves 5.38
P = P_old + 0.5 * (V_old + V) * dt;

end


% %wb, fb, pos(:,i-1), vel(:,i-1), C(:,:,i-1), dti
% function [r_eb_e,v_eb_e,C_b_e] = imuMechECEF(omega_ib_b, f_ib_b, old_r_eb_e, old_v_eb_e, old_C_b_e, tor_i)
% % parameters
% omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)
% 
% % Begins
% 
% % ATTITUDE UPDATE
% % From (2.145) determine the Earth rotation over the update interval
% % C_Earth = C_e_i' * old_C_e_i
% alpha_ie = omega_ie * tor_i;
% C_Earth = [cos(alpha_ie), sin(alpha_ie), 0;...
%           -sin(alpha_ie), cos(alpha_ie), 0;...
%                        0,             0,  1];
%                        
% % Calculate attitude increment, magnitude, and skew-symmetric matrix
% alpha_ib_b = omega_ib_b * tor_i;
% mag_alpha = sqrt(alpha_ib_b' * alpha_ib_b);
% Alpha_ib_b = utils.vec2skew(alpha_ib_b);  
% 
% % Obtain coordinate transformation matrix from the new attitude w.r.t. an
% % inertial frame to the old using Rodrigues' formula, (5.73)
% if mag_alpha>1.E-8
%     C_new_old = eye(3) + sin(mag_alpha) / mag_alpha * Alpha_ib_b +...
%         (1 - cos(mag_alpha)) / mag_alpha^2 * Alpha_ib_b * Alpha_ib_b;
% else
%     C_new_old = eye(3) + Alpha_ib_b;
% end %if mag_alpha    
% 
% % Update attitude using (5.75)
% C_b_e = C_Earth * old_C_b_e * C_new_old;
%     
% % SPECIFIC FORCE FRAME TRANSFORMATION
% % Calculate the average body-to-ECEF-frame coordinate transformation
% % matrix over the update interval using (5.84) and (5.85)
% if mag_alpha>1.E-8
%     ave_C_b_e = old_C_b_e * (eye(3) + (1 - cos(mag_alpha)) / mag_alpha^2 ...
%         * Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha^2 ...
%         * Alpha_ib_b * Alpha_ib_b) - 0.5 * utils.vec2skew([0;0;alpha_ie])...
%         * old_C_b_e;
% else
%      ave_C_b_e = old_C_b_e - 0.5 * utils.vec2skew([0;0;alpha_ie]) *...
%          old_C_b_e;
% end %if mag_alpha     
% 
% % Transform specific force to ECEF-frame resolving axes using (5.85)
% f_ib_e = ave_C_b_e * f_ib_b;
%     
% % UPDATE VELOCITY
% % From (5.36),
% v_eb_e = old_v_eb_e + tor_i * (f_ib_e + utils.gravityECEF(old_r_eb_e) -...
%     2 * utils.vec2skew([0;0;omega_ie]) * old_v_eb_e);
% 
% % UPDATE CARTESIAN POSITION
% % From (5.38),
% r_eb_e = old_r_eb_e + (v_eb_e + old_v_eb_e) * 0.5 * tor_i; 
% 
% end

