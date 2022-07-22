function r = MD_multrnd(p,N)
% _
% Random Numbers from the Multinomial Distribution
% FORMAT r = MD_multrnd(p,N)
% 
%     p - a  1 x K vector with multinomial probabilities
%     N - an integer specifying the number of samples
% 
%     r - an N x 1 vector with multinomial random numbers
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 10:35 (V2)
%  Last edit: 16/10/2018, 10:35 (V2)


% get dimensionality
K = numel(p);

% calculate cumulative probabilities
c = zeros(1,K+1);
for j = 1:K
    c(j+1) = sum(p(1:j));
end;

% generate multinomial random numbers
r = zeros(N,1);
for i = 1:N
    r(i) = sum(rand > c);
end;