function euler = rot2euler(C)
% ROT2EULER converts rotation matrix (nav-to-body) to euler angels
%
% Input:
%   C       3x3 rotaion matrix
%
% Output:
%   euler   3x1 euler angles (roll, pitch, yaw) [rad]
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>

% ZYX not XYZ (Groves 2.23)
phi   =  atan(C(3,2) ./ C(3,3)); % roll
theta = -asin(C(3,1));           % pitch
psi   =  atan2(C(2,1), C(1,1));  % yaw

euler = [phi; theta; psi];

end
