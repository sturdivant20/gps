function euler = quat2euler(q)

% Input:
%   q       4x1 quaternion
%   
% Output:
%   euler   3x1 euler angles (roll, pithc, yaw) [rad]
%
% References:
%   The effects of movement speeds and magnetic disturbance on inerti (2017)
%     -  Howard Chen
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


w = q(1); 
x = q(2); 
y = q(3); 
z = q(4);

% (Chen 2.42)
phi   = atan2(2*(w*x + y*z), (w^2 - x^2 - y^2 + z^2));
theta = asin(-2*(-w*y + x*z));
psi   = atan2(2*(w*z + x*y), (w^2 + x^2 - y^2 - z^2));

euler = [phi; theta; psi];

end
