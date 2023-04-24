function C = quat2rot(q)
% QUAT2ROT transform quaternion into rotation matrix (body-to-nav frame)
%
% Input:
%   q   4x1 quaternion
%
% Output:
%   C   3x3 rotaion matrix
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

% (Chen 2.33)
C(1,1) = w*w + x*x - y*y - z*z;
C(1,2) = 2*(x*y - w*z);
C(1,3) = 2*(w*y + x*z);
C(2,1) = 2*(w*z + x*y);
C(2,2) = w*w - x*x + y*y - z*z;
C(2,3) = 2*(y*z - w*x);
C(3,1) = 2*(x*z - w*y);
C(3,2) = 2*(y*z + w*x);
C(3,3) = w*w - x*x - y*y + z*z;

end
