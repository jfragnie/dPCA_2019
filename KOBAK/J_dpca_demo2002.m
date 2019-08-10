%%%% This section creates toy data.
%
% It should be replaced by actual experimental data. The data should be
% joined in three arrays of the following sizes (for the Romo-like task):
%
% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5) ->returns the sample mean of X, treating NaNs as missing values
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum) 
% % %       -> firingRatesAv.. = mean*size... 

%clear all

% N = 100;%100;   % number of neurons
% T = 20;%20;     % number of time points
% S = 7;%7;       % number of stimuli
% D = 2;          % number of decisions
% E = 20;%20;     % maximal number of trial repetitions
N=28;   %nb of electrodes
T=59;   %nb of time points
E=45;   %nb of trials by rep
D=3;    %nb of brain areas

% 
% 
% % generating three components
%time = (1:T) / 10;
% component{1} = bsxfun(@times, ones(1,S,D,T,E), shiftdim(sin(time), -2)) * 20;
% component{2} = bsxfun(@times, ones(1,S,D,T,E), (1:S)-ceil(S/2)) .* bsxfun(@times, ones(1,S,D,T,E), shiftdim(time, -2));
% component{3} = bsxfun(@times, ones(1,S,D,T,E), shiftdim([-1 1],-1)) .* ...
%     bsxfun(@times, ones(1,S,D,T,E), shiftdim([time(1:end/2)*0 time(1:end/2)], -2)) * 6;
% 
% % mixing them randomly along non-orthogonal axes
% V = randn(N, length(component));
% V(:,2) = V(:,2) + V(:,1)/3;
% V(:,3) = V(:,3) + V(:,2)/4;
% V = bsxfun(@times, V, 1./sqrt(sum(V.^2)));
% for c = 1:length(component)
%     components(:,c) = component{c}(:,:);
% end
% firingRates = V*components';
% %clear component components 
% 
% % adding noise and normalizing
% firingRates = firingRates + randn(size(firingRates));
% firingRates = firingRates - min(firingRates(:)) + 10;
% firingRates = reshape(firingRates, [N S D T E]);
% 
% % setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
%                                  % recordings (they imply the same number 
%                                  % of trials for each neuron)
if ~ifSimultaneousRecording
    trialNum = randi([5 E], [N S D]);
% else  && -> trialNum is created later in "our data" section!
%     %trialNum = repmat(randi([5 E], [1 S D]), [N 1 1]);
end
% 
% for n = 1:N
%     for s = 1:S
%         for d = 1:D
%             firingRates(n,s,d,:,trialNum(n,s,d)+1:E) = nan;
%         end
%     end
% end
% 
% % computing PSTHs
% firingRatesAverage = nanmean(firingRates, 5);

%% our data!
%%load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190212_SmallSphere_200_FiringRate_corr');
%%FR1 = FR_allByEpoch;  ->see J_dpca_demo1202. 
%%clear FR_allByEpoch;
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190220_SmallSphere_200_FiringRate_corr');
FR2 = FR_allByEpoch; 
clear FR_allByEpoch;
load('C:/Users/jules/Downloads/Projet_Bach/HD_data/Rey/Rey2036ElectrodeInfo.mat');
data=H_elec(1:128);
[data_mod,Type] = processData(data);

%% 20.02.19
[indDel2,FR2_red, averageFR] = findSmallIndex(FR2);
[AVG1,AVG2,AVG3] = brainAreas(Type,averageFR);
[firingRates,AVG1,AVG2,AVG3] = dimRed(AVG1,AVG2,AVG3,indDel2,FR2_red);
%trialNum:
[n1,n2,~,~] = size(firingRates); %if 4D
trialNum = size(firingRates,4).*ones(n1,n2);
firingRatesAverage = bsxfun(@times, mean(firingRates,4), size(firingRates,4)./trialNum);

%% Define parameter grouping
% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; here we ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

%does not work anymore with the change I've made! 
% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}; 
% margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - Brain Areas
%    2 - time
%    [1 2] - Brain Areas/time interaction
% They could be grouped as follows: 
combinedParams = {{1, [1 2]},{2}};
margNames = {'Brain areas', 'Condition-independent'};
margColours = [23 100 171; 187 20 25]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
%timeEvents = time(round(length(time)/2));

%not sure about it 
%time = [1:size(firingRates,2)];
time = [1:size(firingRates,3)];
%timeEvents = 1000;
timeEvents = time(round(length(time)/2)); %rounded to 30.

% check consistency between trialNum and firingRates
% for n = 1:size(firingRates,1)
%     for s = 1:size(firingRates,2)
%         for d = 1:size(firingRates,3)
%             assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
%         end
%     end
% end

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ'); %also produces the "economy size"
                          %decomposition. If X is m-by-n with m >= n, then it is
                          %equivalent to svd(X,0). For m < n, only the first m columns 
                          %of V are computed and S is m-by-m.
W = W(:,1:20);  %reduced to 20 components

% minimal plotting
%FIRST PLOT (pca)
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
%SECOND PLOT (still PCA)
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                      ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);
    %get rid of those params:
    %'time', time,                        ...
    %'timeEvents', timeEvents,               ...

%% Step 2: PCA in each marginalization separately

%THIRD PLOT
dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

%fOURTH PLOT
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);
%%
%close all;
%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');
    %'simultaneous', ifSimultaneousRecording, ...
Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: decoding
%does not work with my data 
% decodingClasses = {[(1:D)' (1:D)'], repmat([1:2], [D 1]), [], [(1:D)' (D+(1:D))']};
% %close all;
% accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'numRep', 1, ...        % increase to 100
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
% 
% accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'numRep', 5, ...        % increase to 100
%     'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
% 
% componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16,                ...
%     'componentsSignif', componentsSignif);



%%
%close all;