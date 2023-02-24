function [x,b,P,DOP,i] = gnssPos(x_sv, psr, psr_var, x_hat)
% GNSSPSUEDO Sattelite psuedorangeing function using Iterative Least Squares
%
% Inputs:
%   x_sv    Nx3 ECEF satellite coordinates [m]
%   psr     Nx1 Measured psuedorange [m]
%   psr_var 1x1 or Nx1 Measurement standard deviation
%   x_hat   (optional) initial condition
%
% Outputs:
%   x       3x1 position estimate [m]
%   b       1x1 clock bias [m]
%   P       3x3 estimate covariance
%   i       number of iterations
%
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%

% get number of measurements
N = size(x_sv, 1);

% determine if weighted or unweighted least squares
if length(psr_var) == 1
    W = eye(N);
    wgt = false;
else
    W = inv(diag(psr_var));
    wgt = true;
end

% initial conditions
if nargin < 4
    x_hat = [0;0;0;0];
end
error = Inf;
i = 0;

% least squares iterations
while (error > min(psr_var.^2)*1e-3) && (i < 10)
% while error > 1e-6
    i = i + 1;
    u = x_sv - x_hat(1:3)';
    r = sqrt(sum(u.^2,2));
    uv = u./r;

    H = [-uv, ones(N,1)];       % geometry matrix [-ux/r, -uy/r, -uz/r, 1]
    y = psr - (r + x_hat(4));   % meas. vector [rho - (r + b)]
    dx = inv(H'*W*H)*H'*W*y;

    x_hat = x_hat + dx;
    error = norm(dx);
end
x = x_hat(1:3);
b = x_hat(4);

% calculate covariance
if wgt
    P = inv(H'*W*H);
else
    P = psr_var^2 .* inv(H'*H);
end
DOP = inv(H'*H);

end
