function [obj G Gb] = gaussian(X, y, beta, beta0, dextra)
N = size(X, 1);
obj = norm(X*beta+beta0 - y)^2/N;
G = 2*X'*(X*beta+beta0 - y)/N;
Gb = 2*ones(1, N)*(X*beta+beta0 - y)/N;
end