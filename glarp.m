function sol = glarp(Series, model, init, par, opt, index)
% This function fits the GLARP model to the time series provided in Series
% INPUTS:
%           Series:     A NxT matrix
%           Model:      Specifies (1) the name of the model: model.name and
%                       extra parameters of the model if needed
%           init:       Initialization of the evolution matrices
%           par:        Parameters of the fit, such as lambda and number of
%                       lags used.
%           opt:        Optimization options
%           index:      The indeces that are used for learning: index
%           should not be < par.lags

[N, T] = size(Series);
M = length(index);
P = par.lags;
% opt.verbose = 0;    % To avoid the details of inner iterations
obj = 0;
sol.A = cell(P, 1);
for i = 1:P; sol.A{i} = zeros(N); end
sol.b = zeros(N, 1);
if opt.verboseOut; fprintf('Number of solved cases #: %5d', 0); end
% Build the design matrix for prediction of 'index' samples
% We can make this part parallel
for i = 1:N     % Solving for the i^th time series
    data.X = zeros(M, P*N);
    data.y = Series(i, index)';
    inN = 1:N;
    for j = 1:P
        data.X(:, inN) = Series(:, (index - j))';
        inN = inN + N;
    end
    
    % Initialization
    init2.beta = zeros(P*N, 1);
    for j = 1:P
        init2.beta(N*(j-1)+1:N*j) = init.A{j}(i, :);
    end
    init2.b = init.b(i);    % Will not be used
    
    % Running the APG
    solTemp = APG(model, data, init2, par, opt);
    
    obj = obj + solTemp.obj;
    
    % Converting the results back
    for j = 1:P
        sol.A{j}(i, :) = solTemp.beta( N*(j-1)+1:N*j )';
    end
    sol.b(i) = solTemp.b;
    
    if opt.verboseOut
        fprintf('%c%c%c%c%c%c', 8,8,8,8,8,8);
        fprintf('%5d ', i);
    end
end
if opt.verboseOut; fprintf('\n'); end

% The model selection criteria %% Should be moved to an upper level function
nnZ = 0;
for i = 1:P
    nnZ = nnZ + sum(sum(abs(sol.A{i}) > par.th));
end
sol.aic = obj + 2*nnZ;      % think about the value of '2' here.
sol.bic = obj + nnZ*log(M*N);