function [obj G Gb] = gumbel(X, y, beta, beta0, extra)
sigma = extra.sigma;
d = size(X, 2);
temp = X*beta+beta0;
obj = mean( (y-temp)/sigma + exp(-(y-temp)/sigma) );
G = mean((X/sigma).*repmat(-1 +exp(-(y-temp)/sigma), 1, d), 1)';
Gb = mean(-1+ exp(-(y-temp)/sigma))/sigma;