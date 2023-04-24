function C = euler2rot(euler)
% EULER2ROT converts euler angles to rotation matrix (body-to-nav)
%
% Input:
%   euler   3x1 euler angles (roll, pitch, yaw) [rad]
%
% Output:
%   C       3x3 rotaion matrix
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


phi = euler(1); theta = euler(2); psi = euler(3);

% (Groves 2.22)
C1 = [1, 0       ,  0       ; ...
      0, cos(phi), -sin(phi); ...
      0, sin(phi),  cos(phi)]; 
 
C2 = [ cos(theta), 0, sin(theta); ...
       0         , 1, 0         ; ...
      -sin(theta), 0, cos(theta)];

C3 = [cos(psi), -sin(psi), 0; ...
      sin(psi),  cos(psi), 0; ...
      0       ,  0       , 1];

% ZYX not XYZ
C = C3 * C2 * C1;

end
