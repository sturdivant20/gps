function g = gravityECEF(r)
% GRAVITY_NED calculates gravity vector in the navigation frame
%
% Input:
%   r   3x1 position vector (ECEF)
%
% Output:
%   g   3xN gravity vector in ecef (ECEF) frame
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2012)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>

Ro = 6378137.0000;      % equatorial radius
w_ie = 7.292115e-5;     % Earth rotation rate
mu = 3.986004418e14;    % WGS84 Earth gravitational constant
J2 = 1.082627E-3;       % WGS84 Earth second gravitational constant

% disance from center of earth
R = sqrt(r' * r);

% gravitational accel (Groves 2.142)
z = 5 * (r(3) / R)^2;
y = -mu / R^3 * (r + 1.5*J2*(Ro / R)^2 * [(1 - z) * r(1); ...
                                          (1 - z) * r(2); ...
                                          (3 - z) * r(3)]);

% centripetal accel (Groves 2.133)
g(1:2,1) = y(1:2) + w_ie^2 * r(1:2);
g(3,1) = y(3);


% r_eb_e = r;
% 
% %Parameters
% R_0 = 6378137; %WGS84 Equatorial radius in meters
% mu = 3.986004418E14; %WGS84 Earth gravitational constant (m^3 s^-2)
% J_2 = 1.082627E-3; %WGS84 Earth's second gravitational constant
% omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)
% 
% % Begins
% 
% % Calculate distance from center of the Earth
% mag_r = sqrt(r_eb_e' * r_eb_e);
% 
% % If the input position is 0,0,0, produce a dummy output
% if mag_r==0
%     g = [0;0;0];
%     
% % Calculate gravitational acceleration using (2.142)
% else
%     z_scale = 5 * (r_eb_e(3) / mag_r)^2;
%     gamma = -mu / mag_r^3 *(r_eb_e + 1.5 * J_2 * (R_0 / mag_r)^2 *...
%         [(1 - z_scale) * r_eb_e(1); (1 - z_scale) * r_eb_e(2);...
%         (3 - z_scale) * r_eb_e(3)]);
% 
%     % Add centripetal acceleration using (2.133)
%     g(1:2,1) = gamma(1:2) + omega_ie^2 * r_eb_e(1:2);
%     g(3) = gamma(3);
%     
% end % if

end


