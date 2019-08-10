%% This section creates toy data.
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
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

clear all

%the nomenclature is chosen to be consistent with Kobak's code
N=28;   %nb of electrodes
T=59;   %nb of time points
S=2;    %nb of days
%E      %nb of trials by rep [will be defined later]
D=3;    %nb of brain areas

% setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)

% for n = 1:N
%     for s = 1:S
%         for d = 1:D
%             firingRates(n,s,d,:,trialNum(n,s,d)+1:E) = nan;
%         end
%     end
% end

% computing PSTHs
% firingRatesAverage = nanmean(firingRates, 5);
%% our data!

%data of the electrodes:
load('C:/Users/jules/Downloads/Projet_Bach/HD_data/Rey/Rey2036ElectrodeInfo.mat');
data=H_elec(1:128);
%ordering the electrodes by type.
[data_mod,Type] = processData(data);
clear H_elec data;
FR_struct = {};
E_v = [];  %finding the minimal value for the trials 

%% 12.02.19
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190212_SmallSphere_200_FiringRate_corr');
FR1 = FR_allByEpoch; 
clear FR_allByEpoch;

[~,FR1_red, averageFR1] = findSmallIndex(FR1);
[AVG1_1,AVG2_1,AVG3_1] = brainAreas(Type,averageFR1);
FR_struct{1} = dimRedallDays(AVG1_1,AVG2_1,AVG3_1,FR1_red);
E_v = [E_v length(FR_struct{1,1}(1,1,1,:))];

%% 20.02.19
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190220_SmallSphere_200_FiringRate_corr');
FR2 = FR_allByEpoch; 
clear FR_allByEpoch;

[~,FR2_red, averageFR2] = findSmallIndex(FR2);
[AVG1_2,AVG2_2,AVG3_2] = brainAreas(Type,averageFR2);
FR_struct{2} = dimRedallDays(AVG1_2,AVG2_2,AVG3_2,FR2_red);
E_v = [E_v length(FR_struct{1,2}(1,1,1,:))];

%% number of trials and 5D full matrix
E = min(E_v);
firingRates = zeros(28,3,length(FR_struct),59,E);
for i=length(FR_struct)
   if (length(FR_struct{1,i}(1,1,1,:)) > E) 
      Start=E+1; End=length(FR_struct{1,i}(1,1,1,:));
       FR_struct{1,i}(:,:,:,Start:1:End) = [];   
   end
   %creation of the 5d matrix
   firingRates(:,:,i,:,:) = FR_struct{1,i};
end
[n1,n2,n3,~,~] = size(firingRates); %if 4D, should work with firingRatesX
trialNum = E.*ones(n1,n2,n3);
firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum);

%% clean
clear averageFR1 averageFR2 AVG1_1 AVG1_2 AVG2_1 AVG2_2 AVG3_1 AVG3_2;
clear E_v n1 n2 n3 FR_struct FR1 FR2 Type; 
%% Define parameter grouping
% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - brain areas 
%    2 - days
%    3 - time
% There are three pairwise interactions:
%    [1 3] - brain areas/time interaction
%    [2 3] - days/time interaction
%    [1 2] - brain areas/days interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Brain areas', 'Days', 'Condition-independent', 'B areas/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
time = [0:0.02:1.16];
timeEvents = time(round(length(time)/2));


%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


%% Step 2: PCA in each marginalization separately

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

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);


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
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

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

% decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
% 
% accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'numRep', 5, ...        % increase to 100
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