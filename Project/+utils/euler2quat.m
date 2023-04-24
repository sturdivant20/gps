function q = euler2quat(euler)
% EULER2QUAT converts euler angles to rotation matrix (body-to-nav)
%
% Input:
%   euler   3x1 euler angles (roll, pitch, yaw) [rad]
%
% Output:
%   q       4x1 quaternion
%
% References:
%   The effects of movement speeds and magnetic disturbance on inerti (2017)
%     -  Howard Chen
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


% rearrange to YPR
euler = [euler(3), euler(2), euler(1)];

c = cos( euler/2 );
s = sin( euler/2 );

% ZYX rotation sequence (Chen 2.41)
q = [c(:,1)*c(:,2)*c(:,3) + s(:,1)*s(:,2)*s(:,3); ...
     c(:,1)*c(:,2)*s(:,3) - s(:,1)*s(:,2)*c(:,3); ...
     c(:,1)*s(:,2)*c(:,3) + s(:,1)*c(:,2)*s(:,3); ...
     s(:,1)*c(:,2)*c(:,3) - c(:,1)*s(:,2)*s(:,3)];

end