function p = MD_unipdf(x, a, b)
% _
% Probability Density Function for a Continuous Uniform Distribution
% FORMAT p = MD_unipdf(x, a, b)
% 
%     x - an N x 1 vector of values (domain: a <= x_i <= b)
%     a - a scalar, the lower end of the uniform distribution
%     b - a scalar, the upper end of the uniform distribution
% 
%     p - an N x 1 vector of probability densities
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/08/2018, 13:50 (V1)
%  Last edit: 07/08/2018, 13:50 (V1)


% calculate probability density
p = 1/(b-a) * ones(size(x));