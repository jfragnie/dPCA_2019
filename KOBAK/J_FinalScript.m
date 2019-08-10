%%GENERAL REMARKS
%Elec1 corresponds to the M1 brain area.
%Elec2 corresponds to the PMv (ventral premotor cortex).
%Elec3 corresponds to the PMd (dorsal premotor cortex).

%% laoding our data
% informations about electrodes and the output matrix that will be sorted.

% loading the data of the electrodes
load('C:/Users/jules/Downloads/Projet_Bach/HD_data/Rey/Rey2036ElectrodeInfo.mat');

% loading the matrix with our data.
% 3D matrix with 128x36000x63 doubles, corresponding to 
% firing rate channels x time x trials (i.e movement, the single grasps)
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190212_SmallSphere_200_FiringRate_corr');
smallSphere12_02 = FR_allByEpoch;
clear FR_allByEpoch
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190220_SmallSphere_200_FiringRate_corr');
smallSphere20_02 = FR_allByEpoch;
clear FR_allByEpoch;



%% Processing and sorting of electrodes data
%Substracting the useless 16 arrays of the stuct. Rey2036ElectrodeInfo.
data = H_elec(1:128);
clear H_elec;
[data_mod,Type] = processData(data);

%% Changing the bin size to 20ms! [XXss = data!] 
% group the elements 600 by 600 in the time direction.
% size became dim(data)= [128 60 63] or 59 as they did..

%1st Idea
%for i=1:60
%   XXss(:,i,:) = mean(SmallSphere(:,1+(i-1)*600:i*600,:),2)*10;  %% Unit conversion: to have herz:
%end                                                              %  We have to divise by /0.1 ms.
% mm=mean(mean(XXss,3),2); %for SS 72 val 
% delInd = find(mm<1);

%2nd Idea, use a Gaussian Kernel to smoothen the data
%12.02.2019
[indDel_1,smallSphere1_red, averageFR1] = findSmallIndex(smallSphere12_02);
%20.02.2019
[indDel_2,smallSphere2_red, averageFR2] = findSmallIndex(smallSphere20_02);
clear smallSphere1 smallSphere2;


%% Processing of the problem data
% computing the average of the single grasps (the 3rd dim of our data) over the
% number of trials which is stored in the third dimensions of our 3D matrix.
% recall: nanmean(X,DIM) takes the mean along the dimension DIM of X
%AVGtot = nanmean(FR_allByEpoch,3);  %%with long time vector
%AVGred = nanmean(XXss,3);             %%short concatenated time vector
%AVGred_1d = nanmean(AVGred,2);        %  to be discarded probably
averageAvgFR1 = nanmean(averageFR1,2);
averageAvgFR2 = nanmean(averageFR2,2);
[data_mod,Type] = processData(data);
%12.02.19
[AVG1_1,AVG2_1,AVG3_1] = brainAreas(Type,averageFR1);
%20.02.19
[AVG1_2,AVG2_2,AVG3_2] = brainAreas(Type,averageFR2);

%% PCA
%Now we will first do a PCA for the AVG1,2 &3 and try to compare them.
%numbers 1-3 represents the electrodes to which the average was taken
% /!\ precise to which brain areas they come from when you know it...
% The three first PCs explain sufficient variance so we can discard the
% other
%% DATA 12.02.19
[COEFF1_1,~,~,~, EXPLAINED1_1] = pca(AVG1_1,'NumComponents', 3);
[COEFF2_1,~,~,~, EXPLAINED2_1] = pca(AVG2_1,'NumComponents', 3);
[COEFF3_1,~,~,~, EXPLAINED3_1] = pca(AVG3_1,'NumComponents', 3);
EXPLAINED3_1= [EXPLAINED3_1; zeros(16,1)];
explained_var = [EXPLAINED1_1 EXPLAINED2_1 EXPLAINED3_1];
%Perhaps SCORE will be useful
clear EXPLAINED1_1 EXPLAINED2_1 EXPLAINED3_1;
% We can see that according to the matlab build-in function pca, the explained 
% variances for the three electrodes 1-3 is over 98% taking only the three
% principal components. 
% We can assume then that taking only three PCs is sufficient to explain
% our data with sufficiently high degree of certainty.

%% DATA 20.02.19
[COEFF1_2,~,~,~, EXPLAINED1_2] = pca(AVG1_2,'NumComponents', 3);
[COEFF2_2,~,~,~, EXPLAINED2_2] = pca(AVG2_2,'NumComponents', 3);
[COEFF3_2,~,~,~, EXPLAINED3_2] = pca(AVG3_2,'NumComponents', 3);
EXPLAINED3_2= [EXPLAINED3_2; zeros(16,1)];
explained_var = [EXPLAINED1_2 EXPLAINED2_2 EXPLAINED3_2];
%Perhaps SCORE will be useful
clear EXPLAINED1_2 EXPLAINED2_2 EXPLAINED3_2;



%% Comparison of brain areas
% We will do pair-wise comparison here following the method proposed 
% by Björck and Golub. See equation (1) and (2) in Gallego's paper.
%% DATA 12.02.19
% Comparison of elec1-elec2 area
[~,S1_2before,~] = svd(COEFF1_1'*COEFF2_1);

% Comparison of elec1-elec3 area
[~,S1_3before,~] = svd(COEFF1_1'*COEFF3_1);

% Comparison of elec2-elec3 area
[~,S2_3before,~] = svd(COEFF2_1'*COEFF3_1);
% I don't think that U, V provide informational data for further analysis.
% 1st and 2nd angles are really close to zero.
% 1st angle: cos(theta) > 0.998 for three interactions.
% 2nd angle: cos(theta) > 0.90.
% Interaction 1-3 has the smallest angle: cos(theta) > 0.90 for three
% angles. 
clear COEFF1_1 COEFF2_1 COEFF3_1;

%% DATA 20.02.19
% Comparison of elec1-elec2 area
[~,S1_2after,~] = svd(COEFF1_2'*COEFF2_2);

% Comparison of elec1-elec3 area
[~,S1_3after,~] = svd(COEFF1_2'*COEFF3_2);

% Comparison of elec2-elec3 area
[~,S2_3after,~] = svd(COEFF2_2'*COEFF3_2);
% I don't think that U, V provide informational data for further analysis.
% 1st and 2nd angles are really close to zero.
% 1st angle: cos(theta) > 0.998 for three interactions.
% 2nd angle: cos(theta) > 0.90.
% Interaction 1-3 has the smallest angle: cos(theta) > 0.90 for three
% angles. 
clear COEFF1_2 COEFF2_2 COEFF3_2;

%% QUESTION [for section above; PCA]
% Other way to do it but average it over time I think. is it more correct?
%
% X = bsxfun(@minus, AVG1, mean(AVG1,2));
% 
% [W,~,~] = svd(X, 'econ'); %also produces the "economy size"
%                           %decomposition. If X is m-by-n with m >= n, then it is
%                           %equivalent to svd(X,0). For m < n, only the first m columns 
%                           %of V are computed and S is m-by-m.
% W = W(:,1:20);  %reduced to 20 components


%% dPCA.
% linear dimensionality reduction method, dPCA. This approach identifies a
% single neural manifold for all the data.
%12.02.19 data
[firingRates1,~,~,~,AVG1_1,AVG2_1,AVG3_1] = dimRed(AVG1_1,AVG2_1,AVG3_1,indDel_1,smallSphere1_red);
%20.02.19
[firingRates2,~,~,~,AVG1_2,AVG2_2,AVG3_2] = dimRed(AVG1_2,AVG2_2,AVG3_2,indDel_2,smallSphere2_red);
clear indDel_1 indDel_2;
%[F,D]  = myDPCA(FR_allByEpoch);  %%Should it give information about the marginalization?


