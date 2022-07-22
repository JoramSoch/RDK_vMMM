function [y, b] = ME_preproc(y, b, c)
% _
% Data Pre-Processing for a von Mises Mixture Model
% FORMAT [y, b] = ME_preproc(y, b)
% 
%     y - an N x 1 vector of angles (domain: -pi <= y_i <= +pi)
%     b - a scalar, the systematic error (usually mildly positive)
%     c - a logical, indicating mean-centering (true by default)
% 
%     y - an N x 1 vector of pre-processed angles
%     b - a scalar, the error determined from data
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/10/2018, 12:00 (V2)
%  Last edit: 08/11/2018, 00:35 (V3)


% transform row to column vector, if necessary
if size(y,1) == 1, y = y'; end;

% remove dysfunctional NaN values, if necessary
y(isnan(y)) = [];

% transform degree to radians, if necessary
if max(abs(y)) > (2*pi), y = y * (pi/180); end;

% move to different interval, if necessary
if all(y>0), y = y - pi; end;

% specify systematic error, if necessary
if nargin < 2 || isempty(b), b = mean(y); end;

% specify mean-centering, if necessary
if nargin < 3 || isempty(c), c = true; end;

% mean-center by cyclic subtraction
if c
    y = y - b;
    y(y<-pi) = y(y<-pi) + (2*pi);
    y(y>+pi) = y(y>+pi) - (2*pi);    
end;
