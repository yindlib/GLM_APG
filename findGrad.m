function [obj G Gb] = findGrad(X, y, beta, beta0)
N = size(X, 1);
obj = norm(X*beta+beta0 - y)^2/N;
G = 2*X'*(X*beta+beta0 - y)/N;
Gb = 2*ones(1, N)*(X*beta+beta0 - y)/N;