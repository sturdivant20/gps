function W_en_n = rate_transport(v, r)
% RATE_TRANSPORT transport rate in nav-frame
%
% Input:
%   v       3x1 velocity vector [m/s]
%   r       3x1 position vector [rad, rad, m]
%
% Output:
%   W_en_n  3x3 transport rate matrix (skew symmetric) [rad/s]
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%       Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


vn  = v(1);
ve  = v(2);
L   = r(1);
hgt = r(3);

% meridian and transverse radii of curvature
[Rn, Re] = utils.radii(L);

% angular velocity transport rate (Groves 5.44)
w_en_n = [ ve / (Re + hgt); ...   
          -vn / (Rn + hgt); ...
          -ve * tan(L) / (Re + hgt)];

W_en_n = utils.vec2skew(w_en_n);

end
