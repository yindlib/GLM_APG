function cverr = glarpCV(Series, model, init, par, opt, nCV)
% Cross validation for GLARP

T = size(Series, 2);
P = par.lags;

% Assigning the folds
index = randperm(T-P)+P;
folds = cell(1, nCV);
foldLength = floor((T-P)/nCV);
for i = 1:nCV-1; folds{i} = index(foldLength*(i-1)+1:foldLength*i); end
folds{nCV} = index(foldLength*(nCV-1)+1:end);

% Cross-validation
err = zeros(1, nCV);
optCV = opt;
optCV.verboseOut = 0;   % Silencing the GLARP
for i = 1:nCV
    % Setup the index
    index = [];
    for j = 1:nCV
        if j ~= i
            index = [index folds{j}];
        end
    end
    % Train
    model.sol = glarp(Series, model, init, par, optCV, index);
    % Get the error on the validation fold
    err(i) = pred(Series, model, folds{i});
    err(i) = err(i) * length(folds{i}); % Since the error is normalized in Pred
end
cverr = sum(err);