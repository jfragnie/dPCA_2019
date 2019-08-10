function [firingRates,AVG1tot,AVG2tot,AVG3tot,AVG1,AVG2,AVG3] = dimRed(AVG1,AVG2,AVG3,indDel,FR_allByEpoch_red)
%% General
% From the brain areas obtain in [..]= brainAreas(..), this method will proccess 
%the brain areas obtained in order to have the 3 areas with the same number
% of electrodes. For this we will select data with highest firing rates averaged over time and 
% trials.
% Then we will compute a 4D matrix with those reduced data.
% The matrix is of the form: Electrodes_selected x Time_points x Trials x Brain_areas
% with dimension 28 x 3 x 59 x 63.

%%  Dimensionality reduction to the same size.
%indices to be discarded for each brain areas.
indDel1=find(indDel<=48);  ind = indDel1; indDel1 = indDel(indDel1);
indDel3=find(indDel>96);  ind = [ind; indDel3]; indDel3 = indDel(indDel3);
indDel2 = setdiff(1:length(indDel), ind);  indDel2 = indDel(indDel2);

%indices removed
AVG1(indDel1,:) = [];      %element kept were chosen having the biggest avg firing rate
AVG2(indDel2-48,:) = [];   %substract to have the indices corresponding to the indices of AVG2 (relative information)
AVG3(indDel3-96,:) = [];
AVG1tot=AVG1; 
AVG2tot=AVG2; 
AVG3tot=AVG3;
FR_allByEpoch_red(indDel,:,:) = []; %has to be done here bc indicies will be shift otherwise.

%removing values from AVG1, AVG2 and AVG3 so all three matrices have the same size.                          
AVG1mean=nanmean(AVG1,2);   AVG2mean=nanmean(AVG2,2);
[~, ind1] = sort(AVG1mean, 'descend');    [~, ind2] = sort(AVG2mean, 'descend'); 
siz=size(AVG3,1);
AVG1=AVG1(ind1(1:siz),:);    AVG2=AVG2(ind2(1:siz),:);
ind2=ind2+48; %we need absolute information
newIndDel = [ind1(siz+1:end); ind2(siz+1:end)];  %deletion array with value with smallest firing rate average.
FR_allByEpoch_red(newIndDel,:,:) = [];
%clear AVG1mean AVG2mean ind ind1 ind2 indDel_1 indDel_2 indDel_3 newIndDel indDel i FR_allByEpoch;

%% Create the 4D matrix 
M1 = FR_allByEpoch_red(1:siz,:,:);
PMv = FR_allByEpoch_red(siz+1:2*siz,:,:);
PMd = FR_allByEpoch_red(2*siz+1:end,:,:);
firingRates=zeros(siz,3,59,size(FR_allByEpoch_red,3));
firingRates(:,1,:,:)=M1;
firingRates(:,2,:,:)=PMv;
firingRates(:,3,:,:)=PMd;


end

