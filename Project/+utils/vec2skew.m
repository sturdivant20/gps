function Y = vec2skew(x)
% SKEW skew symmetric form of vector
%
% Input:
%   x   3x1 vector
%
% Output:
%   Y   3x3 skew symmetric matrix
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


% (Groves 2.50)
Y = [    0, -x(3),  x(2); ...
      x(3),     0, -x(1); ...
     -x(2),  x(1),     0];

end