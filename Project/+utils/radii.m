function [Rn, Re] = radii(L)
% RADII_OF_CURVATURE calculates meridian and transverse radii of curvature
%
% Input:
%   L   latitude (radians)
%
% Output:
%   Rn  meridian radius of curvature [m]
%   Re  normal (transverse) radius of curvature [m]
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


Ro = 6378137.0000;      % equatorial radius
Rp = 6356752.3142;      % polar radius

e2 = 1 - (Rp^2/Ro^2);   % eccentricity squared (Groves 2.92)


den = 1 - e2 .* (sin(L)).^2;
Rn = Ro .* (1 - e2) ./ den.^1.5;    % meridian radius (Groves 2.105)
Re = Ro ./ sqrt(den);               % transverse radius (Groves 2.106)

end