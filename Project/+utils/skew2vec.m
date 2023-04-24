function x = skew2vec(Y)
% SKEW2VECTOR transformsskew symmetric form to vector
%
% Input:
%   Y   3x3 skew symmetric matrix
%
% Output:
%   x   3x1 vector
%
% References:
%   Principles of GNSS, Inertial, and Multisensor Integrated Navigation 
%   Systems (2008)  -  Paul D. Groves
%
% Date:     11/2022
% Author:   Daniel Sturdivant <dfs0012@auburn.edu>


% (Groves 2.50)
x = [Y(3,2); Y(1,3); Y(2,1)];

end