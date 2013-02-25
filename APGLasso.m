function [beta b] = APGLasso(X, y, lambda)
ep = 1e-6;
Max_iter = 100;
alpha = 0.8;
delta = 1;

[N d] = size(X);
beta = zeros(d, 1);
objOld = 0;
b = 0;
t = 1;
objs = zeros(1, Max_iter);
for i = 1:Max_iter
    % Perform an iteration to find P_L(Y_k)
    [obj G Gb] = findGrad(X, y, beta, b);
    objs(i) = obj;
    Ybeta = beta - delta * G;
    Ybeta = (abs(Ybeta) > lambda).*(abs(Ybeta)-lambda).*sign(Ybeta);
    Yb = b - delta *Gb;
    
    % Backtracking
    while Obj(X, y, Ybeta, Yb, lambda) > QL(obj, G, Gb, Ybeta, Yb, beta, b, delta, lambda)
        delta = delta * alpha;
        Ybeta = beta - delta * G;
        Yb = b - delta *Gb;
    end
    
    t_new = (1+sqrt(1+4*t^2))/2;
%     t_new = 1;
    beta_new = Ybeta;
    b_new = Yb;
    
    b = b_new + ((t-1)/t_new)*(b_new - b);  
    beta = beta_new + ((t-1)/t_new)*(beta_new - beta);  
    t = t_new;
    
    % Check for termination
    if abs(objOld - obj) < ep
        break
    else
        objOld = obj;
    end
end
plot(objs, 'g')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = QL(obj, G, Gb, Ybeta, Yb, beta, b, delta, lam)
out = obj + (Ybeta - beta)'*G + (Yb - b)*Gb + lam*norm(beta, 1);
out = out + 0.5*norm(Ybeta - beta)^2/delta + 0.5/delta*(b - Yb)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = Obj(X, y, Ybeta, Yb, lam)
obj = findGrad(X, y, Ybeta, Yb);
obj = obj + lam*norm(Ybeta, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obj G Gb] = findGrad(X, y, beta, beta0)
N = size(X, 1);
obj = norm(X*beta+beta0 - y)^2/N;
G = 2*X'*(X*beta+beta0 - y)/N;
Gb = 2*ones(1, N)*(X*beta+beta0 - y)/N;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%