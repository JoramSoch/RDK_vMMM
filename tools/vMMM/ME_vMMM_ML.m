function [r_est, m_est, k_est, MLL] = ME_vMMM_ML(y, mu, kr, m)
% _
% Maximum Likelihood Estimation for a von Mises Mixture Model
% FORMAT [r_est, m_est, k_est, MLL] = ME_vMMM_ML(y, mu, kr, m)
% 
%     y  - an N x 1 vector of angles (domain: -pi <= y_i <= +pi)
%     mu - a  1 x K vector of means (default: mu_1 = 0, mu_2 = pi)
%     kr - a  1 x 2 vector of precision range (default: kr = [1 21])
%     m  - a string indicating the model to be used for estimation
% 
%     r_est - a 1 x K vector of estimated frequencies
%     m_est - a 1 x K vector of estimated locations
%     k_est - a 1 x K vector of estimated precisions
%     MLL   - a scalar, the maximum log-likelihood
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/08/2018, 14:15 (V1)
%  Last edit: 08/11/2018, 00:50 (V3)


% set mu and kappa
if nargin < 2 || isempty(mu), mu = [0 pi]; end;
if nargin < 3 || isempty(kr), kr = [1 21]; end;
if nargin < 4 || isempty(m) , m  = 'm11';  end;

% m00: the null model ("only guessing")
if strcmp(m,'m00') || strcmp(m,'m00b');
    r_est = [0 0 1];
    m_est = [0 pi];
    k_est = [0 0];
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m01: extended model ("guessing and ROOD"), w/o bias [old V7a]
if strcmp(m,'m01')
    % set initial values
    r  = [0, 0.5, 0.5];
    k  = [0, (kr(2)+1)/2];
    dr = 0.25;
    dk = (kr(2)-kr(1))/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r2 = [(r(2)-2*dr):dr:(r(2)+2*dr)];
        k2 = [(k(2)-2*dk):dk:(k(2)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k(1));
        Ik2 = besseli(0,k2);
        % log-likelihood
        LL = zeros(numel(r2),numel(k2));
        for i2 = 1:numel(r2)
            rs = [r(1), r2(i2), (1-r2(i2))];
            for j2 = 1:numel(k2)
                ks = [k(1), k2(j2)];
                LL(i2,j2) = ME_vMMM_LL(y, rs, mu, ks, [Ik1, Ik2(j2)]);
            end;
        end;
        % adjust parameters
        [MLL, ind] = max(LL(:));
        [r2i, k2i] = ind2sub(size(LL), ind);
        r = [r(1), r2(r2i), (1-r2(r2i))];
        k = [k(1), k2(k2i)];
        % adjust step sizes
        dr = dr/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = mu;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m01: extended model ("guessing and ROOD"), with bias [new V2]
if strcmp(m,'m01b')
    % set initial values
    r  = [0, 0.5, 0.5];
    b  = [0];
    k  = [0, mean(kr)];
    dr = r(2)/2;
    db = pi/8;
    dk = range(kr)/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r2 = [(r(2)-2*dr):dr:(r(2)+2*dr)];
        mb = [(b-2*db):db:(b+2*db)];
        k2 = [(k(2)-2*dk):dk:(k(2)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k(1));
        Ik2 = besseli(0,k2);
        % log-likelihood
        LL = zeros(numel(r2),numel(mb),numel(k2));
        for i2 = 1:numel(r2)
            rs = [r(1), r2(i2), (1-r2(i2))];
            for ij = 1:numel(mb)
                ms = [mu(1), mu(2)+mb(ij)];
                for j2 = 1:numel(k2)
                    ks = [k(1), k2(j2)];
                    LL(i2,ij,j2) = ME_vMMM_LL(y, rs, ms, ks, [Ik1, Ik2(j2)]);
                end;
            end;
        end;
        % adjust parameters
        [MLL, ind]      = max(LL(:));
        [r2i, mbi, k2i] = ind2sub(size(LL), ind);
        r = [r(1), r2(r2i), (1-r2(r2i))];
        b =  mb(mbi);
        m = [mu(1), mu(2)+mb(mbi)];
        k = [k(1), k2(k2i)];
        % adjust step sizes
        dr = dr/2;
        db = db/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = m;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m10: extended model ("guessing and detection"), w/o bias [old V7a]
if strcmp(m,'m10')
    % set initial values
    r  = [0.5, 0, 0.5];
    k  = [(kr(2)+1)/2, 0];
    dr = 0.25;
    dk = (kr(2)-kr(1))/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r1 = [(r(1)-2*dr):dr:(r(1)+2*dr)];
        k1 = [(k(1)-2*dk):dk:(k(1)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k1);
        Ik2 = besseli(0,k(2));
        % log-likelihood
        LL = zeros(numel(r1),numel(k1));
        for i1 = 1:numel(r1)
            rs = [r1(i1), r(2), (1-r1(i1))];
            for j1 = 1:numel(k1)
                ks = [k1(j1), k(2)];
                LL(i1,j1) = ME_vMMM_LL(y, rs, mu, ks, [Ik1(j1), Ik2]);
            end;
        end;
        % adjust parameters
        [MLL, ind] = max(LL(:));
        [r1i, k1i] = ind2sub(size(LL), ind);
        r = [r1(r1i), r(2), (1-r1(r1i))];
        k = [k1(k1i), k(2)];
        % adjust step sizes
        dr = dr/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = mu;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m10: extended model ("guessing and detection"), with bias [new V2]
if strcmp(m,'m10b')
    % set initial values
    r  = [0.5, 0, 0.5];
    b  = [0];
    k  = [mean(kr), 0];
    dr = r(1)/2;
    db = pi/8;
    dk = range(kr)/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r1 = [(r(1)-2*dr):dr:(r(1)+2*dr)];
        mb = [(b-2*db):db:(b+2*db)];
        k1 = [(k(1)-2*dk):dk:(k(1)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k1);
        Ik2 = besseli(0,k(2));
        % log-likelihood
        LL = zeros(numel(r1),numel(mb),numel(k1));
        for i1 = 1:numel(r1)
            rs = [r1(i1), r(2), (1-r1(i1))];
            for ij = 1:numel(mb)
                ms = [mu(1)+mb(ij), mu(2)];
                for j1 = 1:numel(k1)
                    ks = [k1(j1), k(2)];
                    LL(i1,ij,j1) = ME_vMMM_LL(y, rs, ms, ks, [Ik1(j1), Ik2]);
                end;
            end;
        end;
        % adjust parameters
        [MLL, ind]      = max(LL(:));
        [r1i, mbi, k1i] = ind2sub(size(LL), ind);
        r = [r1(r1i), r(2), (1-r1(r1i))];
        b =  mb(mbi);
        m = [mu(1)+mb(mbi), mu(2)];
        k = [k1(k1i), k(2)];
        % adjust step sizes
        dr = dr/2;
        db = db/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = m;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m11: full model ("guessing + ROOD + detection"), w/o bias [old V7a]
if strcmp(m,'m11')
    % set initial values
    r  = [0.5, 0.5, 0];
    k  = [(kr(2)+1)/2, (kr(2)+1)/2];
    dr = 0.25;
    dk = (kr(2)-kr(1))/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r1 = [(r(1)-2*dr):dr:(r(1)+2*dr)];
        r2 = [(r(2)-2*dr):dr:(r(2)+2*dr)];
        k1 = [(k(1)-2*dk):dk:(k(1)+2*dk)];
        k2 = [(k(2)-2*dk):dk:(k(2)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k1);
        Ik2 = besseli(0,k2);
        % log-likelihood
        LL = zeros(numel(r1),numel(r2),numel(k1),numel(k2));
        for i1 = 1:numel(r1)
            for i2 = 1:numel(r2)
                if r2(i2) > (1-r1(i1))
                    LL(i1,i2,:,:) = -Inf;
                else
                    rs = [r1(i1), r2(i2), (1-r1(i1)-r2(i2))];
                    for j1 = 1:numel(k1)
                        for j2 = 1:numel(k2)
                            ks = [k1(j1), k2(j2)];
                            LL(i1,i2,j1,j2) = ME_vMMM_LL(y, rs, mu, ks, [Ik1(j1), Ik2(j2)]);
                        end;
                    end;
                end;
            end;
        end;
        % adjust parameters
        [MLL, ind]           = max(LL(:));
        [r1i, r2i, k1i, k2i] = ind2sub(size(LL), ind);
        r = [r1(r1i), r2(r2i), (1-r1(r1i)-r2(r2i))];
        k = [k1(k1i), k2(k2i)];
        % adjust step sizes
        dr = dr/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = mu;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;

% m11: full model ("guessing + ROOD + detection"), with bias [new V2]
if strcmp(m,'m11b')
    % set initial values
    r  = [0.5, 0.5, 0];
    b  = [0];
    k  = [mean(kr), mean(kr)];
    dr = r(1)/2;
    db = pi/8;
    dk = range(kr)/4;
    cc = 1e-4;
    % refined grid search
    while dr > cc
        % parameter grid
        r1 = [(r(1)-2*dr):dr:(r(1)+2*dr)];
        r2 = [(r(2)-2*dr):dr:(r(2)+2*dr)];
        mb = [(b-2*db):db:(b+2*db)];
        k1 = [(k(1)-2*dk):dk:(k(1)+2*dk)];
        k2 = [(k(2)-2*dk):dk:(k(2)+2*dk)];
        % Bessel function
        Ik1 = besseli(0,k1);
        Ik2 = besseli(0,k2);
        % log-likelihood
        LL = zeros(numel(r1),numel(r2),numel(mb),numel(k1),numel(k2));
        for i1 = 1:numel(r1)
            for i2 = 1:numel(r2)
                if r2(i2) > (1-r1(i1))
                    LL(i1,i2,:,:) = -Inf;
                else
                    rs = [r1(i1), r2(i2), (1-r1(i1)-r2(i2))];
                    for ij = 1:numel(mb)
                        ms = mu + mb(ij);
                        for j1 = 1:numel(k1)
                            for j2 = 1:numel(k2)
                                ks = [k1(j1), k2(j2)];
                                LL(i1,i2,ij,j1,j2) = ME_vMMM_LL(y, rs, ms, ks, [Ik1(j1), Ik2(j2)]);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        % adjust parameters
        [MLL, ind]                = max(LL(:));
        [r1i, r2i, mbi, k1i, k2i] = ind2sub(size(LL), ind);
        r = [r1(r1i), r2(r2i), (1-r1(r1i)-r2(r2i))];
        b = mb(mbi);
        m = mu + mb(mbi);
        k = [k1(k1i), k2(k2i)];
        % adjust step sizes
        dr = dr/2;
        db = db/2;
        dk = dk/2;
    end;
    % collect final estimates
    r_est = r;
    m_est = m;
    k_est = k;
    MLL   = ME_vMMM_LL(y, r_est, m_est, k_est);
end;