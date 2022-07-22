function p = ME_vMMM_pdf(y, r, mu, ka, Ik)
% _
% Probability Density Function for a von Mises Mixture Model
% FORMAT p = ME_vMMM_pdf(y, r, mu, ka, Ik)
% 
%     y  - an N x 1 vector of response deviations
%     r  - a  1 x K vector of mixture frequencies
%     mu - a  1 x K vector of means
%     ka - a  1 x K vector of precisions
%     Ik - a  1 x K vector of Bessel values
% 
%     p  - an N x 1 vector of probability densities
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 11:15 (V2)
%  Last edit: 16/10/2018, 11:15 (V2)


% if frequencies don't add up to 1
if abs(sum(r)-1) > 1e-10
    p = zeros(size(y));
else
    % if frequencies are not in [0,1]
    if any(r<0) || any(r>1)
        p = zeros(size(y));
    else
        % if precisions are smaller than 0
        if any (ka<0)
            p = zeros(size(y));
        else
            % otherwise, calculate densities
            if nargin < 5 || isempty(Ik)
                p = r(1) * MD_vmpdf(y, mu(1), ka(1)) + r(2) * MD_vmpdf(y, mu(2), ka(2)) + r(3) * MD_unipdf(y, -pi, +pi);
            else
                p = r(1) * MD_vmpdf(y, mu(1), ka(1), Ik(1)) + r(2) * MD_vmpdf(y, mu(2), ka(2), Ik(2)) + r(3) * MD_unipdf(y, -pi, +pi);
            end;
        end;
    end;
end;