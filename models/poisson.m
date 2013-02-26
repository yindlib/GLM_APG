function [obj G Gb] = poisson(X, y, beta, beta0, dextra)
d = size(X, 2);
temp = X*beta+beta0;
obj = -mean(y.*temp - exp(temp));
G = -mean(X.*repmat(y - exp(temp), 1, d), 1)';
Gb = -mean(y-exp(temp));