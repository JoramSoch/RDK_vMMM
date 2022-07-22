function r = MD_vmrnd(mu, ka, N)
% _
% Random Numbers from a von Mises Distribution
% FORMAT r = MD_vmrnd(mu, ka, N)
% 
%     mu - the mean of the von Mises distribution
%     ka - the precision of the von Mises distribution
%     N  - the number of samples to be drawn
% 
%     r  - an N x 1 vector of von Mises random numbers
% 
% FORMAT r = MD_vmrnd(mu, ka, N) returns N random numbers r from the
% von Mises distribution with mean mu and precision ka.
% 
% This function uses an envelope-rejection method based on a wrapped Cauchy
% distribution to draw random variates from an arbitrary von Mises
% distribution, first proposed in [1].
% 
% References:
% [1] Best DJ, Fisher NI (1979): "Efficient Simulation of the von Mises
%     Distribution". Journal of the Royal Statistical Society. Series C
%     (Applied Statistics), vol. 28, no. 2, pp. 152-157.
% 
% Author: Dylan Muir <muir@hifo.uzh.ch>, June 19th, 2012.
% Source: https://de.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin
% 
% Editor: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 10:20 (V2)
%  Last edit: 16/10/2018, 10:20 (V2)


% uniform special case
if ka == 0
    r = MD_unirnd(a, b, N);

% non-uniform general case
else % ka > 0
    
    % pre-allocate results
    r = NaN(N,1);
    Z = NaN(N,1);
    F = NaN(N,1);
    C = NaN(N,1);
    
    % pre-compute values
    tau = sqrt(4 * ka^2 + 1) + 1;
    rho = (tau - sqrt(2 * tau)) / (2 * ka);
    R   = (1 + rho^2) / (2 * rho);
    
    % draw random variates
    tbDraw = true(N,1);
    tbAccept = ~tbDraw;
    while any(tbDraw)
        
        % draw partial variates, estimate wrapped Cauchy envelope
        nd = sum(tbDraw);
        Z(tbDraw) = cos(pi * rand(nd,1));
        F(tbDraw) = (1 + R * Z(tbDraw)) ./ (R + Z(tbDraw));
        C(tbDraw) = ka * (R - F(tbDraw));
        
        % filter partial variates
        R2 = rand(nd,1);
        tbAccept(tbDraw) = ((C(tbDraw) .* (2 - C(tbDraw))) - R2) > 0;
        tbRecheck = ~tbAccept(tbDraw);
        if any(tbRecheck)
            tbAccept(tbDraw & ~tbAccept) = (log(C(tbDraw & ~tbAccept) ./ R2(tbRecheck)) + 1 - C(tbDraw & ~tbAccept)) >= 0;
        end;
        
        % construct final variates
        na = sum(tbDraw & tbAccept);
        r(tbDraw & tbAccept) = sign(rand(na,1)-0.5) .* acos(F(tbDraw & tbAccept));
        tbDraw(tbAccept) = false;
        
    end;
    
    % shift to mean parameter
    r = r + mu;
    r(r<-pi) = r(r<-pi) + (2*pi);
    r(r>+pi) = r(r>+pi) - (2*pi);
    
end;