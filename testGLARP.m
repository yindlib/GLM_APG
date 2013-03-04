clc
addpath('./models')
T = 10000;
N = 40;
model.dextra.sigma = 1;
model.dextra.nu = 1;
%% Solver settings
opt.tol = 1e-4;
opt.Max_iter = 500;
opt.alpha = 0.8;
opt.delta = 1;
opt.verbose = 1;    % Don't show inner iterations
opt.verboseOut = 1; % Show the details of the outer iterations
opt.nob = 1;

init.b = zeros(N, 1);
init.A = cell(1, 1);
init.A{1} = zeros(N);

par.lambda = 0.00001;
par.lags = 1;
%% Gaussian case
A1 = [0.9, 0, 0, 0; 1, 0.9, 0, 0; 1, 0, 0.9, 0; 1, 0, 0, 0.9];
A = kron(eye(N/4), A1);
Series = zeros(N, T);
Series(:, 1) = randn(N, 1);
for t = 2:T
    Series(:, t) = A*Series(:, t-1) + 0.1*randn(N, 1);
end

index = 2:T;

model.fname = 'gaussian';
sol = glarp(Series, model, init, par, opt, index);

% AGPLassoL(Series, 'smth')