function r = MD_unirnd(a, b, N)
% _
% Random Numbers from a Continuous Uniform Distribution
% FORMAT r = MD_unirnd(a, b, N)
% 
%     a - the lower end of the uniform distribution
%     b - the upper end of the uniform distribution
%     N - the number of samples to be drawn
% 
%     r  - an N x 1 vector of uniform random numbers
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 09:55 (V2)
%  Last edit: 16/10/2018, 09:55 (V2)


% generate random numbers
r = (rand(N,1) * (b-a)) + a;