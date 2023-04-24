function [P, V, C] = imuMechNAV(w_ib_b, f_ib_b, P_old, V_old, C_old, dt)
% NAV_MECH IMU body-to-nav frame mechanization
%
% Inputs:
%   w_ib_b  3x1 IMU angular velcoity reading in body frame [rad/s]
%   f_ib_b  3x1 IMU specific force readings in body frame [m/s^2]
%   P_old   3x1 previous position vector LLA [rad, rad, m]
%   V_old   3x1 previous velocity vector in nav-frame [m/s]
%   C_old   3x3 previous body-to-nav rotation matrix
%   dt      incremental time [s]
%
% Outputs:
%   P       3x1 updated position [rad, rad, m]
%   V       3x1 updated velcoity [m/s]
%   C       3x3 updated rotation matrix
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%       Systems (2013)  -  Paul D. Groves
%


%% VELOCITY UPDATE

% old rotation rates
W_en_n = utils.rateTransport(V_old, P_old);
W_ie_n = utils.rateEarth(P_old);

% propagated angular rates
Alpha = utils.vec2skew(w_ib_b * dt);

% specific force transformation (Groves 5.83/5.84/5.86)
C_avg = C_old * expm(Alpha./2) - 0.5*(W_ie_n + W_en_n)*C_old*dt;
f_ib_n = C_avg * f_ib_b;

% velocity update (Groves 5.54)
g = utils.gravityNED(P_old);
V = V_old + (f_ib_n + g - (W_en_n + 2*W_ie_n)*V_old)*dt;


%% POSITION UPDATE

% old meridian radius of curvature
[Rn_o, Re_o] = utils.radii(P_old(1));

% (Groves 5.56)
hgt_n = P_old(3) - 0.5*(V_old(3) + V(3))*dt;
lat_n = P_old(1) + 0.5*((V_old(1) / (Rn_o + P_old(3))) + ...
                        (V(1)     / (Rn_o + hgt_n  )))*dt;

% new meridian radius of curvature
[~, Re_n] = utils.radii(lat_n);

lon_n = P_old(2) + 0.5*((V_old(2) / (cos(P_old(1))*(Re_o + P_old(3)))) + ...
                        (V(2)     / (cos(lat_n)   *(Re_n + hgt_n   ))))*dt;

% position update
P = [lat_n; lon_n; hgt_n];


%% ATTITUDE UPDATE

% new rotation rates
W_en_n_new = utils.rateTransport(V, P);

% Groves 5.73 (Rodrigues' Formula)
C_bar = expm(Alpha);

% update attitude (Groves 5.76)
C = C_old*C_bar - (W_ie_n + 0.5*W_en_n + 0.5*W_en_n_new)*C_old*dt;

end

% % function [P, V, C] = imuMechNAV(w_ib_b, f_ib_b, P_old, V_old, C_old, dt)
% function [P, v_eb_n, C_b_n] = imuMechNAV(...
%     omega_ib_b, f_ib_b, P_old, old_v_eb_n, old_C_b_n, tor_i)
% %Nav_equations_NED - Runs precision local-navigation-frame inertial 
% %navigation equations (Note: only the attitude update and specific force
% %frame transformation phases are precise.)
% %
% % Software for use with "Principles of GNSS, Inertial, and Multisensor
% % Integrated Navigation Systems," Second Edition.
% %
% % This function created 1/4/2012 by Paul Groves
% %
% % Inputs:
% %   tor_i         time interval between epochs (s)
% %   old_L_b       previous latitude (rad)
% %   old_lambda_b  previous longitude (rad)
% %   old_h_b       previous height (m)
% %   old_C_b_n     previous body-to-NED coordinate transformation matrix
% %   old_v_eb_n    previous velocity of body frame w.r.t. ECEF frame, resolved
% %                 along north, east, and down (m/s)
% %   f_ib_b        specific force of body frame w.r.t. ECEF frame, resolved
% %                 along body-frame axes, averaged over time interval (m/s^2)
% %   omega_ib_b    angular rate of body frame w.r.t. ECEF frame, resolved
% %                 about body-frame axes, averaged over time interval (rad/s)
% % Outputs:
% %   L_b           latitude (rad)
% %   lambda_b      longitude (rad)
% %   h_b           height (m)
% %   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
% %                 north, east, and down (m/s)
% %   C_b_n         body-to-NED coordinate transformation matrix
% 
% % Copyright 2012, Paul Groves
% % License: BSD; see license.txt for details
% 
% old_L_b = P_old(1);
% old_lambda_b = P_old(2);
% old_h_b = P_old(3);
% 
% % parameters
% omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)
% 
% % Begins
% 
% % PRELIMINARIES
% % Calculate attitude increment, magnitude, and skew-symmetric matrix
% alpha_ib_b = omega_ib_b * tor_i;
% mag_alpha = sqrt(alpha_ib_b' * alpha_ib_b);
% Alpha_ib_b = utils.vec2skew(alpha_ib_b);  
% 
% % From (2.123) , determine the angular rate of the ECEF frame
% % w.r.t the ECI frame, resolved about NED
% omega_ie_n = omega_ie * [cos(old_L_b); 0; - sin(old_L_b)];
%     
% % From (5.44), determine the angular rate of the NED frame
% % w.r.t the ECEF frame, resolved about NED
% [old_R_N,old_R_E] = utils.radii(old_L_b);
% old_omega_en_n = [old_v_eb_n(2) / (old_R_E + old_h_b);...
%     -old_v_eb_n(1) / (old_R_N + old_h_b);...
%     -old_v_eb_n(2) * tan(old_L_b) / (old_R_E + old_h_b)];
%     
% % SPECIFIC FORCE FRAME TRANSFORMATION
% % Calculate the average body-to-ECEF-frame coordinate transformation
% % matrix over the update interval using (5.84) and (5.86)
% if mag_alpha>1.E-8
%     ave_C_b_n = old_C_b_n * (eye(3) + (1 - cos(mag_alpha)) / mag_alpha^2 ...
%         * Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha^2 ...
%         * Alpha_ib_b * Alpha_ib_b) -...
%         0.5 * utils.vec2skew(old_omega_en_n + omega_ie_n) * old_C_b_n;
% else
%      ave_C_b_n = old_C_b_n -...
%          0.5 * utils.vec2skew(old_omega_en_n + omega_ie_n) * old_C_b_n;
% end %if mag_alpha     
% 
% % Transform specific force to ECEF-frame resolving axes using (5.86)
% f_ib_n = ave_C_b_n * f_ib_b;
%     
% % UPDATE VELOCITY
% % From (5.54),
% v_eb_n = old_v_eb_n + tor_i * (f_ib_n + utils.gravityNED(P_old) -...
%     utils.vec2skew(old_omega_en_n + 2 * omega_ie_n) * old_v_eb_n);
% 
% % UPDATE CURVILINEAR POSITION
% % Update height using (5.56)
% h_b = old_h_b - 0.5 * tor_i * (old_v_eb_n(3) + v_eb_n(3));
% 
% % Update latitude using (5.56)
% L_b = old_L_b + 0.5 * tor_i * (old_v_eb_n(1) / (old_R_N + old_h_b) +...
%     v_eb_n(1) / (old_R_N + h_b));
% 
% % Calculate meridian and transverse radii of curvature
% [R_N,R_E]= utils.radii(L_b);
% 
% % Update longitude using (5.56)
% lambda_b = old_lambda_b + 0.5 * tor_i * (old_v_eb_n(2) / ((old_R_E +...
%     old_h_b) * cos(old_L_b)) + v_eb_n(2) / ((R_E + h_b) * cos(L_b))); 
% 
% P = [L_b; lambda_b; h_b];
% 
% % ATTITUDE UPDATE
% % From (5.44), determine the angular rate of the NED frame
% % w.r.t the ECEF frame, resolved about NED
% omega_en_n = [v_eb_n(2) / (R_E + h_b);...
%     -v_eb_n(1) / (R_N + h_b);...
%     -v_eb_n(2) * tan(L_b) / (R_E + h_b)];
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
% % Update attitude using (5.77)
% C_b_n = (eye(3) - utils.vec2skew(omega_ie_n + 0.5 * omega_en_n + 0.5 *...
%     old_omega_en_n)  * tor_i) * old_C_b_n * C_new_old;

