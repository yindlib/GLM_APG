clc
N = 500;
d = 50;
X = rand(N, d);
beta = [1 1 1 zeros(1, d-3)]';
beta0 = 0;
y = X*beta +beta0 +0.2*randn(N, 1);

data.X = X;
data.y = y;

opt.tol = 5e-3;
opt.Max_iter = 300;
opt.alpha = 0.8;
opt.delta = 1;
opt.verbose = 2;

init.b = 0;
init.beta = zeros(d, 1);

par.lambda = 0.000001;
fname = 'gaussian';
sol = APG(fname, data, init, par, opt)