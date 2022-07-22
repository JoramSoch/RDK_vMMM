function [tt, thABc, pp, thABj] = ME_vMMM_trl(y, r_est, k_est)
% _
% Trial-Wise Analysis for a von Mises Mixture Model
% FORMAT [tt, thABc, pp, thABj] = ME_vMMM_trl(y, r_est, k_est)
% 
%     y     - an N x 1 vector of response deviations
%     r_est - a  1 x K vector of estimated frequencies
%     k_est - a  1 x K vector of estimated precisions
% 
%     tt    - an N x 1 vector of discrete trial types (trial classification)
%     thABc - a  1 x 2 vector of trial type bounds (conditional probabilities)
%     pp    - an N x 3 matrix of posterior probabilities (trial inference)
%     thABj - a  1 x 2 vector of trial type bounds (joint probabilities)
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 11:30 (V2)
%  Last edit: 16/10/2018, 11:30 (V2)


% get number of trials
N = numel(y);

% calculate conditional bounds
thA =  acos((1/k_est(1)) * log(besseli(0,k_est(1))));
thB = -acos((1/k_est(2)) * log(besseli(0,k_est(2)))) + pi;
thABc = [thA, thB];

% determine trial types
tt = zeros(N,1);
tt(y >= -thA & y <= +thA) = 1;
tt(y <= -thB | y >= +thB) = 2;
tt(tt == 0) = 3;

% calculate posterior probalities
pp = zeros(N,3);
pp(:,1) = (r_est(1) * MD_vmpdf(y, 0, k_est(1))) ./ ME_vMMM_pdf(y, r_est, [0 pi], k_est);
pp(:,2) = (r_est(2) * MD_vmpdf(y, pi, k_est(2))) ./ ME_vMMM_pdf(y, r_est, [0 pi], k_est);
pp(:,3) = (r_est(3) / (2*pi)) ./ ME_vMMM_pdf(y, r_est, [0 pi], k_est);

% calculate joint bounds
thA =  acos((1/k_est(1)) * log((r_est(3)/r_est(1)) * besseli(0,k_est(1))));
thB = -acos((1/k_est(2)) * log((r_est(3)/r_est(2)) * besseli(0,k_est(2)))) + pi;
thABj = [thA, thB];