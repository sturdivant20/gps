function [C,R] = ned2ecef(lla)
% NED2ECEF creates rotation matrix from NED to ECEF
%
% Input:
%   lla     3x1 global coordinates [rad, rad, m]
%
% Output:
%   R       3x3 rotaion ned-to-ecef matrix
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>

cos_lat = cos(lla(1));
sin_lat = sin(lla(1));
cos_long = cos(lla(2));
sin_long = sin(lla(2));

% R = [-sin(lla(1,1))*cos(lla(2,1)), -sin(lla(2,1)), -cos(lla(1,1))*cos(lla(2,1)); ...
%      -sin(lla(1,1))*sin(lla(2,1)),  cos(lla(2,1)), -cos(lla(1,1))*sin(lla(2,1)); ...
%                     cos(lla(1,1)),              0,               -sin(lla(1,1))];

% ecef2ned
R = [-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
               -sin_long,            cos_long,        0;...
     -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];

% ned2ecef
C = R';
end

