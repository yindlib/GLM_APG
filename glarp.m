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

% Build the design matrix for prediction of 'index' samples
data.X = zeros(M * N, N*(P*N + 1) );   % N time series and N*L+1 samples for each.
data.y = zeros(M * N, 1);

inT = 1:M;
for i = 1:N
    data.y(inT) = Series(i, index);
    inN = (N*(i-1)+1):N*i;
    for j = 1:P
        data.X(inT, inN) = Series(:, (index - j))';
        inN = inN + N^2;
    end
    data.X(inT, end-(N-i)) = 1;
    inT = inT + M;
end

% Initialization
init2.beta = zeros(N*(P*N + 1), 1);
for j = 1:P
    init2.beta(N^2*(j-1)+1:N^2*j) = reshape(init.A{j}', 1, N^2);
end
init2.beta(end-N+1:end) = init.b';
init2.b = 0;    % Will not be used
opt.nob = 1;

% Running the APG
solTemp = APG(model, data, init2, par, opt);

% Converting the results back
sol.A = cell(P, 1);
for j = 1:P
    sol.A{j} = reshape(solTemp.beta( N^2*(j-1)+1:N^2*j ), N, N)';
end
sol.b = solTemp.beta(end-N+1:end)';
sol.obj = solTemp.obj;