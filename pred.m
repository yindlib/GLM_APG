function [err normerr] = pred(Series, model, index)
% This function fits the GLARP model to the time series provided in Series
% INPUTS:
%           Series:     A NxT matrix
%           Model:      Specifies (1) the name of the model: model.name (2)
%                       extra parameters of the model if needed (3) model
%                       parameters
%           index:      The indeces that are used for learning: index
%           should not be < par.lags
M = length(index);

% First find the eta in GLM for each prediction point
eta = repmat(model.sol.b, 1, M);
for t = 1:length(index)
    for ll = 1:length(model.sol.A)
        eta(:, t) = eta(:, t) + model.sol.A{ll}*Series(:, index(t)-ll);
    end
end

% Compute the MSE
tsmean = mean(mean( abs(Series(:, index)) ));   % For normalized prediction error
switch model.fname
    case 'gaussian'
        err = norm(eta-Series(:, index), 'fro')/numel(eta);
    case 'poisson'
        err = norm(exp(eta)-Series(:, index), 'fro')/numel(eta);
    case 'comp' % Veryfy this
        err = norm(exp(eta+ 1./model.dextra.nu)-Series(:, index), 'fro')/numel(eta);
    case 'gumbel'
        err = norm(eta -psi(1)-Series(:, index), 'fro')/numel(eta);
    otherwise
        err = 0;
        error('Misspecified model.')
end
normerr = err/tsmean;