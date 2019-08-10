function [Xtest, Xtrain] = dpca_getTestTrials(firingRatesPerTrial, numOfTrials, varargin)
%firingRatesPerTrial == firingRates
%optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
%dpca_optimizeLambda(Xfull, Xtrial, numOfTrials, varargin)
%[Xtest, XtrainFull] = dpca_getTestTrials(Xtrial, numOfTrials, 'simultaneous', options.simultaneous);
%[Xtest, XtrainFull] = dpca_getTestTrials(firingRates, trialNum, 'simultaneous', false);


options = struct('simultaneous', false);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = size(firingRatesPerTrial);  %[84    59    63     3]

if ~options.simultaneous
    neuronsConditions = numOfTrials(:);  %with their data: [N S D]
    testTrials = ceil(rand([length(neuronsConditions) 1]) .* neuronsConditions);
    
else
    neuronsConditions = numOfTrials(:);
    neuronsConditions = neuronsConditions(1:size(numOfTrials,1):end); 
    testTrials = ceil(rand([length(neuronsConditions) 1]) .* neuronsConditions);
    testTrials = bsxfun(@times, ones(size(numOfTrials,1),1), testTrials');
    testTrials = testTrials(:);
    
end

ind = reshape(testTrials, size(numOfTrials));  %ind has the size of numTrials : in our case 84x3.
ind = bsxfun(@times, ones(dim(1:end-1)), ind);  %reshape it to a 4D string with dim of firingRatesAVG
%ind = bsxfun(@times, ones(84,59,3), ind);  
%as our data do not have the same shape we have to change the size of ones(dim(1:end-1))

ind = ind(:);

indtest = sub2ind([prod(dim(1:end-1)) dim(end)], (1:prod(dim(1:end-1)))', ind);
%indtest = sub2ind([prod(84,3) 59], (1:prod(84,3))', ind);


Xtest = firingRatesPerTrial(indtest);
Xtest = reshape(Xtest, dim(1:end-1));

% Note to myself (12 Jul 2016): This function used to return only Xtest,
% and all statistics of Xtrain were updated without ever constructing
% Xtrain explicitly. I am now changing this for simplicity (in order to
% implement simultaneously recorded case), but this requries 2x memory.
if nargout > 1
    Xtrain = firingRatesPerTrial;
    Xtrain(indtest) = nan;
end
