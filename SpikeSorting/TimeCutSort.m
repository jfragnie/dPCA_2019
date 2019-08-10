function combinedOverTrials = TimeCutSort(H,unitTimes,NEV)
TimeCut = {};
Length = -0.2*30000:30000;
d1 = size(unitTimes,1);
d2 = size(Length,2);
d3 = length(H.Data.indMovNSP1);
combinedOverTrials = zeros(d1,d2,d3);
%unitTimes = one of the unitTimes matrices!

for i=1:length(H.Data.indMovNSP1)
    if (H.Data.indMovNSP1(i)>NEV.MetaTags.DataDuration && H.Data.indMovNSP1(i)+30000>NEV.MetaTags.DataDuration)
        break
    end
    TimeCut{i} = unitTimes(:,H.Data.indMovNSP1(i)-0.2*30000:H.Data.indMovNSP1(i)+1*30000)/30000;
    
   % combinedOverTrials(:,:,i) = TimeCut{i};
end

end