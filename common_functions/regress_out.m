function Y_res = regress_out(Y, X1)
% Regress X1 out of X2 (both Nx1 or NxK)
% Returns residuals of X2 after removing linear effect of X1

% ensure column vectors
X1 = X1(:);

% design matrix with intercept
X = [ones(size(X1)), X1];

% solve least squares
beta = X \ Y;

% residuals
Y_res = Y - X * beta;
end