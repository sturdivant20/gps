function [xDot, bDot] = gnssVel(x_sv, v_sv, x_user, dopp)
% GNSSPSUEDO calculates user velocity given doppler and satellite
% positions using Least Squares
%
% Inputs:
%   x_sv        Nx3 ECEF satellite coordinates [m]
%   x_user      3x1 ECEF position estimate of the user [m]
%   v_sv        Nx3 ECEF satellite velocity [m]
%   dopp        Nx1 Measured doppler [m]
%   dopp_var    1x1 or Nx1 Measurement standard deviation
%
% Outputs:
%   xDot    3x1 ECEF velocity estimate [m]
%   bDot    1x1 clock drift [m]
%   PDot    3x3 estimate covariance
%
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%

N = size(x_sv, 1);
i = 0;
x = [0;0;0;0];
error = Inf;

% while (error > 1e-6) && (i < 10)
%     i = i+1;

    % generate unit vector and range to satellite
    u = x_sv - x_user';
    r = sqrt(sum(u.^2,2));
    uv = u./r;
    
    % least squares parameters
    H = [-uv, ones(N,1)];
    y = dopp - (sum(uv.*v_sv, 2) + x(4));
    dx = inv(H'*H)*H'*y;

    x = x + dx;
    error = norm(dx);
% end

xDot = x(1:3);
bDot = x(4);

end

