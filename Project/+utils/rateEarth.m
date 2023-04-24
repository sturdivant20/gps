function W_ie_n = rateEarth(L)
% RATE_EARTH rotation rate of earth in nav-frame
%
% Input:
%   L       latitude [rad]
%
% Output:
%   W_ie_n  3x3 earth rate matrix (skew symmetric) [rad/s]
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%       Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>

% angular velocity of earth (Groves 2.123)
w_ie_n = 7.292115e-5 * [cos(L); 0; -sin(L)];
W_ie_n = utils.vec2skew(w_ie_n);

end