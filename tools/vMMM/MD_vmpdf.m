function p = MD_vmpdf(x, mu, ka, Ik)
% _
% Probability Density Function for a von Mises Distribution
% FORMAT p = MD_vmpdf(x, mu, ka)
% 
%     x  - an N x 1 vector of angles (domain: -pi <= x_i <= +pi)
%     mu - a scalar, the mean of the von Mises distribution
%     ka - a scalar, the precision of the von Mises distribution
%     Ik - a scalar, the modified Bessel function of order 0
% 
%     p  - an N x 1 vector of probability densities
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/08/2018, 13:50 (V1)
%  Last edit: 07/08/2018, 13:50 (V1)


% calculate probability density
if nargin < 4 || isempty(Ik)
    p = exp(ka * cos(x-mu)) / (2*pi * besseli(0,ka));
else
    p = exp(ka * cos(x-mu)) / (2*pi * Ik);
end;