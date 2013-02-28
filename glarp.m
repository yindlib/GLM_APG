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
%           index:      The indeces that are used for learning

[N, T] = size(Series);
M = length(index);
P = par.lags;

% Build the design matrix for prediction of 'index' samples
Am = zeros(M * N, N*(P*N + 1) );   % N time series and N*L+1 samples for each.
bm = zeros(M * N, 1);

for j = 1:N     % Add this part to the matrix
    for i = index
        bm(N*(j-1) + i) = Series(j, i);
        Am(N*(j-1) + i, :) = reshape([fliplr(Series(:, (i-P):(i-1))), ones(N, 1)]', 1, N*(P+1));
    end
end

data.X = Am;
data.y = bm;

% Convert the initiaization, too.

solTemp = APG(model.fname, data, init, par, opt, model.dextra);

% Convert the solution