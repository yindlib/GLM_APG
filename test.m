clc
addpath('./models')
N = 500;
d = 50;
X = rand(N, d);
beta = [1 1 1 zeros(1, d-3)]';
beta0 = 0;
model.dextra.sigma = 1;
model.dextra.nu = 1;
%% Solver settings
opt.tol = 1e-4;
opt.Max_iter = 500;
opt.alpha = 0.8;
opt.delta = 1e-2;
opt.verbose = 2;
opt.nob = 1;

init.b = 0;
init.beta = zeros(d, 1);

par.lambda = 0.001;
%% Gaussian case
y = X*beta +beta0 +0.2*randn(N, 1);

data.X = X;
data.y = y;

model.fname = 'gaussian';
sol = APG(model, data, init, par, opt);
%% Poisson case
y = poissrnd(X*beta+beta0);
data.X = X;
data.y = y;

model.fname = 'poisson';
sol = APG(model, data, init, par, opt);
%% COM-Poisson case
fname = 'comp';
init.b = 0.01;
sol = APG(model, data, init, par, opt);
%% Gumbel case
mu = X*beta+beta0;
y = mu - (model.dextra.sigma).*log(-log(rand(size(mu))));
data.X = X;
data.y = y;

model.fname = 'gumbel';
sol = APG(model, data, init, par, opt);
