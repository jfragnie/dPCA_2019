function firingRates = dimRedallDays(AVG1,AVG2,AVG3,FR_allByEpoch_red)
%% General
% From the brain areas obtain in [..]= brainAreas(..), this method will proccess 
%the brain areas obtained in order to have the 3 areas with the same number
% of electrodes. For this we will select data with highest firing rates averaged over time and 
% trials.
% Then we will compute a 4D matrix with those reduced data.
% The matrix is of the form: Electrodes_selected x Time_points x Trials x Brain_areas
% with dimension 28 x 3 x 59 x 63.

%%  Dimensionality reduction to the same size.

M1_flat = nanmean(AVG1,2);
PMv_flat = nanmean(AVG2,2);
PMd_flat = nanmean(AVG3,2);
%take the 28 smallest comp:
[~,  ind1] = sort(M1_flat, 'descend');
[~, ind2] = sort(PMv_flat, 'descend');
[~, ind3] = sort(PMd_flat, 'descend');
ind1red = ind1(1:28); diff1 = setdiff(ind1,ind1red); 
ind2red = ind2(1:28); diff2 = setdiff(ind2,ind2red);
ind3red = ind3(1:28); diff3 = setdiff(ind3,ind3red);
deleteInd = [diff1; diff2+48; diff3+96];

AVG1 = AVG1(ind1,:);
AVG2 = AVG2(ind2,:);
AVG3 = AVG3(ind3,:);
FR_allByEpoch_red(deleteInd,:,:) = []; %has to be done here bc indicies will be shift otherwise.


%% Create the 4D matrix 
M1 = FR_allByEpoch_red(1:28,:,:);
PMv = FR_allByEpoch_red(29:56,:,:);
PMd = FR_allByEpoch_red(57:end,:,:);
firingRates=zeros(28,3,59,size(FR_allByEpoch_red,3));
firingRates(:,1,:,:)=M1;
firingRates(:,2,:,:)=PMv;
firingRates(:,3,:,:)=PMd;


end

