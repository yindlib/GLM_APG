clc
addpath('./models')
N = 500;
d = 50;
X = rand(N, d);
beta = [1 1 1 zeros(1, d-3)]';
beta0 = 0;
dextra.sigma = 1;
dextra.nu = 1;
%% Solver settings
opt.tol = 1e-4;
opt.Max_iter = 500;
opt.alpha = 0.8;
opt.delta = 1;
opt.verbose = 2;

init.b = 0;
init.beta = zeros(d, 1);

par.lambda = 0.001;

%% Gaussian case
y = X*beta +beta0 +0.2*randn(N, 1);

data.X = X;
data.y = y;

fname = 'gaussian';
sol = APG(fname, data, init, par, opt, dextra);

%% Poisson case
y = poissrnd(X*beta+beta0);
data.X = X;
data.y = y;

fname = 'poisson';
sol = APG(fname, data, init, par, opt, dextra);

%% COM-Poisson case
fname = 'comp';
init.b = 0.01;
sol = APG(fname, data, init, par, opt, dextra);
%% Gumbel case
mu = X*beta+beta0;
y = mu - (dextra.sigma).*log(-log(rand(size(mu))));
data.X = X;
data.y = y;

fname = 'gumbel';
sol = APG(fname, data, init, par, opt, dextra);