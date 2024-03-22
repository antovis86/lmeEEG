function [tval, b, se] = lmeEEG_regress(y,X)
% y can be a column [n 1] or matrix [n, nY] now, to allow it to work across 
% all time-points at once. edited by John P. Grogan, 2024
% adapted from:
% REGRESS Multiple linear regression using least squares.
%   References:
%      [1] Chatterjee, S. and A.S. Hadi (1986) "Influential Observations,
%          High Leverage Points, and Outliers in Linear Regression",
%          Statistical Science 1(3):379-416.
%      [2] Draper N. and H. Smith (1981) Applied Regression Analysis, 2nd
%          ed., Wiley.
%   Copyright 1993-2014 The MathWorks, Inc.

nY = size(y,2); % number of columns of Y passed in - e.g. time-samples

[n,ncolX] = size(X);
% Use the rank-revealing QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);
if isempty(R)
    p = 0;
elseif isvector(R)
    p = double(abs(R(1))>0);
else
    p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
end
if p < ncolX
    warning(message('stats:regress:RankDefDesignMat'));
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
% Column of coefficients per column in Y
b = zeros(ncolX,nY);
b(perm,:) = R \ (Q'*y);


% compute SE e tval
RI = R\eye(p);
nu = max(0,n-p);                % Residual degrees of freedom, same for all columns
yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
normr = vecnorm(r);             % get norm per column
if nu ~= 0
    rmse = normr/sqrt(nu);      % Root mean square error.
else
    rmse = NaN;
end

se = zeros(ncolX,nY); % per column
se(perm,:) = rmse .* sqrt(sum(abs(RI).^2,2)); % per column
tval = b./se;


end 
