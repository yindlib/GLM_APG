function [obj G Gb] = comp(X, y, beta, beta0, dextra)
nu = ones(size(y))*dextra.nu;
d = size(X, 2);
temp = exp(X*beta+beta0);
obj = -mean(nu.*y.*log(temp) - funcS(temp, nu));
mup = funcSp(temp, nu);
G = -mean(X.*repmat(nu.*y - mup.*temp, 1, d), 1)';
Gb = -mean(nu.*y - mup.*temp);