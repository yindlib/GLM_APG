clc
N = 500;
d = 50;
X = rand(N, d);
beta = [1 1 1 zeros(1, d-3)]';
beta0 = 0;
y = X*beta +beta0 +0.2*randn(N, 1);

[obj G Gb] = findGrad(X, y, beta, beta0);

beta2= beta;
delta = 1e-7;
i=4;
beta2(i) = beta2(i) + delta;
[obj2] = findGrad(X, y, beta2, beta0);

[obj3] = findGrad(X, y, beta, beta0+delta);

disp((obj2-obj)/delta)
disp(G(i))

disp((obj3-obj)/delta)
disp(Gb)

APGLasso(X, y, 0.000002)