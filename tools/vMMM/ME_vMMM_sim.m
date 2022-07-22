function y = ME_vMMM_sim(r, mu, ka, N)
% _
% Simulate Date from a von Mises Mixture Model
% FORMAT y = ME_vMMM_Sim(r, mu, ka, N)
% 
%     r  - a  1 x K vector of mixture frequencies
%     mu - a  1 x K vector of means
%     ka - a  1 x K vector of precisions
% 
%     y  - an N x 1 vector of response deviations
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 10:40 (V2)
%  Last edit: 16/10/2018, 10:40 (V2)


% set r, mu, kappa, N
if nargin < 1 || isempty(r),  r  = [1/3, 1/3, 1/3]; end;
if nargin < 2 || isempty(mu), mu = [0, pi];         end;
if nargin < 3 || isempty(ka), ka = [1, 1];          end;
if nargin < 4 || isempty(N),  N  = 1000;            end;

% sample trial types
t = MD_multrnd(r, N);

% sample response deviations
y = NaN(N,1);
y(t==1) = MD_vmrnd(mu(1), ka(1), sum(t==1));
y(t==2) = MD_vmrnd(mu(2), ka(2), sum(t==2));
y(t==3) = MD_unirnd(-pi, +pi, sum(t==3));