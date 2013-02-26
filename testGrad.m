% Testing Gradients

clc
N = 500;
d = 50;
X = rand(N, d);
beta = [1 1 1 zeros(1, d-3)]';
beta0 = 0;

extra.sigma = 0.5;
extra.nu = 3;

y = X*beta +beta0 +0.2*randn(N, 1);

[obj1 G Gb] = comp(X, y, beta, beta0, extra);

delta = 1e-8;
[obj2] = comp(X, y, beta, beta0+delta, extra);

beta2 = beta;
beta2(2) = beta(2) + delta;

[obj3] = comp(X, y, beta2, beta0, extra);

disp((obj2-obj1)/delta)
disp(Gb)

disp((obj3-obj1)/delta)
disp(G(2))