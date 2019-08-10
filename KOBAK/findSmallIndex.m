function [Index, FR_allByEpoch_red, averageFR] = findSmallIndex(data)
load('C:/Users/jules/Downloads/Projet_Bach/MyProject/data/20190212_SmallSphere_200_FiringRate_corr');

Fs = 30000;
PreTimeLength=0.2*Fs; %Observing time window for the average of the brain (1 sec) after kuka
PostTimeLength=0.4*Fs; %Observing time window for the average of the brain (1 sec) before kuka  %%NOT the same as in PCA_NSP1 (0.2 & 0.3)
PostTimeLength=1*Fs; %Observing time window for the average of the brain (1 sec) before kuka    %%  Overwriting?



FR_allByEpochNEW = data(:,1:PreTimeLength+PostTimeLength,:);
clear data;
FR_allByEpoch = FR_allByEpochNEW;

H.Fs = 30000;
slidingWindow = 40e-3*H.Fs;  %%20ms
overlapWindow = slidingWindow/2;
win1 = 1:slidingWindow:size(FR_allByEpoch,2);
win1= win1-1; %%changed by me.
win2 = overlapWindow:slidingWindow:size(FR_allByEpoch,2);  %%as before

win = unique([win1, win2]);
for tr = 1:size(FR_allByEpoch,3)
    for wn = 1:length(win)-2
        prov(:,wn) = mean(FR_allByEpoch(:,win(wn)+1:win(wn+2),tr),2);  %thus removed '-1' of win(wn+2)-1 bc line17
                                                                       %+1 added to the first term
                                                                       
    end
    prov(:,wn+1) = mean(FR_allByEpoch(:,win(wn+1):end,tr),2);
    FR_allByEpoch_red(:,:,tr) = prov;
    clear prov
end

averageFR = mean(FR_allByEpoch_red,3)*10;   %then average it over the third dimension
                                            %needed to have the mean in the correct units: /0.1ms
%FR_sum = sum(averageFR,2);  %not needed
FR_mean = mean(averageFR,2);  %needed to have the mean in the correct units: /0.1ms   
%FR_peak = max(averageFR,[],2)/0.1; %not needed

Index = find(FR_mean < 1); 