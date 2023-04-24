function [g, g0] = gravityNED(r)
% GRAVITY_NED calculates gravity vector in the navigation frame
%
% Input:
%   r   3x1 position vector
%
% Output:
%   g   3xN gravity vector in nav (NED) frame
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


L = r(1);
hgt = r(3);

N = max(size(L));
g = zeros(3,N);

Ro = 6378137.0000;      % equatorial radius
Rp = 6356752.3142;      % polar radius
e2 = 1 - (Rp^2/Ro^2);   % eccentricity squared (Groves 2.92)
w_ie = 7.292115e-5;     % Earth rotation rate
mu = 3.986004418e14;    % WGS84 Earth gravitational constant
f = (Ro - Rp) / Ro;     % WGS84 flattening (Groves 2.92)

% somigliana model (Groves 2.134)
sinL2 = sin(L).^2;
g0 = 9.7803253359 .* ((1 + 0.001931853.*sinL2) ./ sqrt(1 - e2.*sinL2));

% north gravity (groves 2.140)
g(1,:) = -8.08E-9 .* hgt .* sin(2 .* L);

% east gravity
g(2,:) = 0;

% down gravity
g(3,:) = g0 .* (1 - (2 ./ Ro) .* (1 + f .* (1 - 2 .* sinL2) + ...
    (w_ie.^2 .* Ro.^2 .* Rp ./ mu)) .* hgt + (3 .* hgt.^2 ./ Ro.^2));

end
