function sol = APG(fname, data, init, par, opt, dextra)
% This function finds the minimum of the function provided in 'fname' via
% Accelerated Proximal Gradient method. It has backtracking and termination
% implemented.
% INPUTS:
%       fname:  The name of the function to be optimized.
%       data:   Usually has two components: data.X for features and data.y
%               for responses.
%       init:   Initial values for the optimization variables
%       par:    Hyper parameters, such as lambda in Lasso
%       opt:    Options for the optimization such as Max_iter, initial
%               delta, backtraking factor alpha, tolerance and verbose
%       dextra: The extra parametes that a distribution may need
delta = opt.delta;
beta = init.beta;
objOld = 0;
b = init.b;
t = 1;
objs = zeros(1, opt.Max_iter);
if opt.verbose; fprintf('Iter #: %5d', 0); end
for i = 1:opt.Max_iter
    % Perform an iteration to find P_L(Y_k)
    [obj G Gb] = feval(fname, data.X, data.y, beta, b, dextra);
    objs(i) = obj;
    Ybeta = beta - delta * G;
    Ybeta = (abs(Ybeta) > par.lambda).*(abs(Ybeta)-par.lambda).*sign(Ybeta);
    Yb = b - delta *Gb;
    
    % Backtracking
    while Obj(fname, data.X, data.y, Ybeta, Yb, par.lambda, dextra) > QL(obj, G, Gb, Ybeta, Yb, beta, b, delta, par.lambda)
        delta = delta * opt.alpha;
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
    if abs(objOld - obj) < opt.tol * obj
        break
    else
        objOld = obj;
    end
    
    if opt.verbose
        fprintf('%c%c%c%c%c%c', 8,8,8,8,8,8);
        fprintf('%5d ', i);
    end
end

sol.b = b;
sol.beta = beta;

if opt.verbose 
    fprintf('\n')
end
if opt.verbose > 1
    plot(objs)
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = QL(obj, G, Gb, Ybeta, Yb, beta, b, delta, lam)
out = obj + (Ybeta - beta)'*G + (Yb - b)*Gb + lam*norm(beta, 1);
out = out + 0.5*norm(Ybeta - beta)^2/delta + 0.5/delta*(b - Yb)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = Obj(fname, X, y, Ybeta, Yb, lam, dextra)
obj = feval(fname, X, y, Ybeta, Yb, dextra);
obj = obj + lam*norm(Ybeta, 1);
end