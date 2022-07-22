function [r_est, m_est, k_est, AICw] = ME_vMMM_MA(y, mu, kr, s)
% _
% (Quasi-Bayesian) Model Averaging across von Mises Mixture Models
% FORMAT [r_est, m_est, k_est, AICw] = ME_vMMM_MA(y, mu, kr, s)
% 
%     y  - an N x 1 vector of angles (domain: -pi <= y_i <= +pi)
%     mu - a  1 x K vector of means (default: mu_1 = 0, mu_2 = pi)
%     kr - a  1 x 2 vector of precision range (default: kr = [1 21])
%     s  - an M x 1 cell array of models to be used for averaging
% 
%     r_est - a  1 x K vector of estimated frequencies
%     m_est - a  1 x K vector of estimated locations
%     k_est - a  1 x K vector of estimated precisions
%     AICw  - an M x 1 vector of Akaike weights
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/04/2020, 15:35 (V5)
%  Last edit: 07/04/2020, 15:35 (V5)


% set mu and kappa
if nargin < 2 || isempty(mu), mu = [0 pi]; end;
if nargin < 3 || isempty(kr), kr = [1 21]; end;

% set model space
if nargin < 4 || isempty(s)
    s = {'m00b', 'm01b', 'm10b', 'm11b'}';
end;

% preallocate parameters
M     = numel(s);
r_est = zeros(M,3);
m_est = zeros(M,2);
k_est = zeros(M,2);
MLL   = zeros(M,1);
p     = zeros(M,1);

% estimate parameters
for i = 1:M                     % MLEs of model parameters
    [r_est(i,:), m_est(i,:), k_est(i,:), MLL(i)] = ME_vMMM_ML(y, mu, kr, s{i});
    switch s{i}                 % number of model parameters,
        case 'm00',  p(i) = 0;  % see "ME_vMMM_ML.m"
        case 'm01',  p(i) = 2;
        case 'm10',  p(i) = 2;
        case 'm11',  p(i) = 4;
        case 'm00b', p(i) = 0;
        case 'm01b', p(i) = 3;
        case 'm10b', p(i) = 3;
        case 'm11b', p(i) = 5;
    end;
end;

% assess model quality
AIC  = -2*MLL + 2*p;
AICw = exp(-1/2*(AIC - min(AIC))) ./ sum( exp(-1/2*(AIC - min(AIC))) );

% average parameters
r_est = sum( r_est .* repmat(AICw, [1 size(r_est,2)]) );
m_est = sum( m_est .* repmat(AICw, [1 size(m_est,2)]) );
k_est = sum( k_est .* repmat(AICw, [1 size(k_est,2)]) );