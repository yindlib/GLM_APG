function sol = glarp_par(Series, model, init, par, opt, index)
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

N = size(Series, 1);
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
parfor i = 1:N     % Solving for the i^th time series
    X = zeros(M, P*N);
    y = Series(i, index)';
    inN = 1:N;
    for j = 1:P
        X(:, inN) = Series(:, (index - j))';
        inN = inN + N;
    end
    
    % Initialization
    beta = zeros(P*N, 1);
    for j = 1:P
        beta(N*(j-1)+1:N*j) = init.A{j}(i, :);
    end
    b = init.b(i);    % Will not be used
    
    % Running the APG
    solTemp = specialAPG(model, X, y, beta, b, par, opt);
    
    obj = obj + solTemp.obj;
    
    % Converting the results back
    for j = 1:P
        A(i, :, j) = solTemp.beta( N*(j-1)+1:N*j )';
    end
    b(i) = solTemp.b;
end
if opt.verboseOut; fprintf('\n'); end

for j = 1:P
    sol.A{j} = squeeze(A(:, :, j));
end
% sol.b = b;

% The model selection criteria %% Should be moved to an upper level function
nnZ = 0;
for i = 1:P
    nnZ = nnZ + sum(sum(abs(sol.A{i}) > par.th));
end
sol.aic = obj + 2*nnZ;      % think about the value of '2' here.
sol.bic = obj + nnZ*log(M*N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxilary function to cope with matlab's parfor
function Sol = specialAPG(model, X, y, beta, b, par, opt)
data.X = X;
data.y = y;
init2.beta = beta;
init2.b = b;
Sol = APG(model, data, init2, par, opt);
end