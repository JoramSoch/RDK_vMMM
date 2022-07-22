function [r_CI, m_CI, k_CI] = ME_vMMM_CI(y, r_est, m_est, k_est, m, alpha)
% _
% Confidence Intervals for a von Mises Mixture Model using Wilks's Theorem
% FORMAT [r_CI, m_CI, k_CI] = ME_vMMM_CI(y, r_est, m_est, k_est, m, alpha)
% 
%     y     - an N x 1 vector of angles (domain: -pi <= y_i <= +pi)
%     r_est - a  1 x K vector of estimated frequencies
%     m_est - a  1 x K vector of estimated locations
%     k_est - a  1 x K vector of estimated precisions
%     m     - a string indicating the model to be used for estimation
%     alpha - the significance level, one minus the confidence level
% 
%     r_CI  - a 2 x K matrix of frequency confidence intervals
%     m_CI  - a 2 x K matrix of location confidence intervals
%     k_CI  - a 2 x K matrix of precision confidence intervals
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 06/11/2018, 15:30 (V3)
%  Last edit: 08/11/2018, 01:20 (V3)


% compute Wilks's log-likelihood threshold
MLL    = ME_vMMM_LL(y, r_est, m_est, k_est);
chi2cv = chi2inv(1-alpha,1);
LL_CI  = MLL-(1/2)*chi2cv;

% preallocate CI for parameters
if strcmp(m,'m00') || strcmp(m,'m00b')
    r_CI = [0, 0, 1;
            0, 0, 1];
    m_CI = [0, pi;
            0, pi];
    k_CI = [0, 0;
            0, 0];
else
    r_CI = zeros(2,numel(r_est));
    m_CI = zeros(2,numel(m_est));
    k_CI = zeros(2,numel(k_est));
end;

% estimate CI for location bias
if isempty(strfind(m,'b')) || strcmp(m,'m00b')
    m_CI = [0, pi;
            0, pi];
end;
if ~isempty(strfind(m,'b')) && ~strcmp(m,'m00b')
    CI_upp = get_bound(0, +pi/8, 'mb');
    CI_low = get_bound(0, -pi/8, 'mb');
    m_CI = [m_est+CI_upp; m_est+CI_low];
end;

% estimate CI for frequency (1)
if strncmp(m,'m10',3) || strncmp(m,'m11',3)
    CI_upp = get_bound(r_est(1), +1e-1, 'r1');
    CI_low = get_bound(r_est(1), -1e-1, 'r1');
    r_CI(:,1) = [CI_upp; CI_low];
end;

% estimate CI for frequency (2)
if strncmp(m,'m01',3) || strncmp(m,'m11',3)
    CI_upp = get_bound(r_est(2), +1e-1, 'r2');
    CI_low = get_bound(r_est(2), -1e-1, 'r2');
    r_CI(:,2) = [CI_upp; CI_low];
end;

% estimate CI for frequency (3)
if ~strncmp(m,'m00',3)
    CI_upp = get_bound(r_est(3), +1e-1, 'r3');
    CI_low = get_bound(r_est(3), -1e-1, 'r3');
    r_CI(:,3) = [CI_upp; CI_low];
end;

% estimate CI for precision (1)
if strncmp(m,'m10',3) || strncmp(m,'m11',3)
    CI_upp = get_bound(k_est(1), +1, 'k1');
    CI_low = get_bound(k_est(1), -1, 'k1');
    k_CI(:,1) = [CI_upp; CI_low];
end;

% estimate CI for precision (2)
if strncmp(m,'m01',3) || strncmp(m,'m11',3)
    CI_upp = get_bound(k_est(2), +1, 'k2');
    CI_low = get_bound(k_est(2), -1, 'k2');
    k_CI(:,2) = [CI_upp; CI_low];
end;

% function: identify CI bound
function CI = get_bound(p, dp, pid)

% get ML estimates
rs = r_est;
ms = m_est;
ks = k_est;
cc = 1e-4;                      % convergence criterion
mi = 1e+4;                      % maximum iterations
ll = true;                      % likelihood larger?
% refined grid search
 i = 0;
while abs(dp) > cc && i <= mi
    % update parameter
    p = p + dp;
    switch pid
        case 'r1', rs = [p, r_est(2)-(p-r_est(1))/2, r_est(3)-(p-r_est(1))/2];
        case 'r2', rs = [r_est(1)-(p-r_est(2))/2, p, r_est(3)-(p-r_est(2))/2];
        case 'r3', rs = [r_est(1)-(p-r_est(3))/2, r_est(2)-(p-r_est(3))/2, p];
        case 'mb', ms =  m_est + p;
        case 'k1', ks = [p, k_est(2)];
        case 'k2', ks = [k_est(1), p];
    end;
    % inspect likelihood
    LL = ME_vMMM_LL(y, rs, ms, ks);
    if ll && LL < LL_CI
        dp = dp*(-1/2);
        ll = false;
    elseif ~ll && LL > LL_CI
        dp = dp*(-1/2);
        ll = true;
    end;
    i = i + 1;
end;
% collect final estimate
CI = p;

end;

end