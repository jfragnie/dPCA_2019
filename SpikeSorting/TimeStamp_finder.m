function H = TimeStamp_finder(analog,kinematic,NEV,indMovAdjBeg)   
    for iA = 1:length(analog.labels)
        if strcmp(analog.labels(iA),'Voltage.BlackRock NSP1')
        %if strcmp(analog.labels(iA),'Voltage.BlackRock Trigger')    
            ind_NSP1 = iA;
        end
    end

    iniNSP1_1000 = find(diff(analog.data(:,ind_NSP1))>2)+1;
    iniNSP1_100 = round(iniNSP1_1000/analog.framerate*kinematic.framerate);

    indMovNSP1 = indMovAdjBeg-iniNSP1_100;
    %indMovNSP1 = indMovAdj-iniNSP1_100; % When using Diana's code
    neg = [];
    neg = find(indMovNSP1<0);
    if ~isempty(neg)
        indMovNSP1(neg) = [];
    end
    H.Data.Fs = NEV.MetaTags.SampleRes;  %30'000Hz.
    H.Data.indMovNSP1 = round(indMovNSP1/kinematic.framerate*30000);

    H.Data.indMovNSP1Time=H.Data.indMovNSP1/30000;
end