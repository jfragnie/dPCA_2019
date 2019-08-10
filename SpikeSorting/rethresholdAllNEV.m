% %% Set up the location of the data
% root_dir = 'W:\Data\Experiments';
% monkey = 'Cersei';
% the_date = '20180613';
%
% data_dir = fullfile (root_dir, monkey, the_date, 'BLACKROCK');

threshold_file = 5; % which file index to use for the thresholding
threshold_factor = -6; % do spike identification as threshold_factor * RMS

%% find all of the available .ns5 files
all_files  = dir(fullfile(data_dir,'*.ns5'));


%% load the file to use for thresholding
%fullfile: build full file name from parts.
fpath = fullfile(all_files(threshold_file).folder,all_files(threshold_file).name);
ns = openNSxCervical(fpath);


%% find threshold for each channel using the reference file
n_chans = size(ns.Data,1);
chan_thresholds = zeros(1,n_chans);

fs = ns.MetaTags.SamplingFreq;
fc = [1000 5000]; % cutoff frequency in hz

for chan = 1:n_chans
    disp(['Channel ' num2str(chan) ' of ' num2str(n_chans)]);
    
    signal_raw = double(ns.ElectrodesInfo(1).MaxAnalogValue)/double(ns.ElectrodesInfo(1).MaxDigiValue) * double(ns.Data(chan,:));
    [b,a] = butter(2,fc/(fs/2));
    signal_filt = filtfilt(b,a,signal_raw);
    
    chan_thresholds(chan) = threshold_factor * rms(signal_filt);
end

%%
% chan_thresholds = [NEV.ElectrodesInfo.LowThreshold];
% chan_thresholds = chan_thresholds(1:130);

%% loop along all files and rethreshold
for iFile = 1:length(all_files)
    disp(['File ' num2str(iFile) ' of ' num2str(length(all_files))]);
    file_name = all_files(iFile).name(1:end-4);
    file_name_save = ['ReTh_' file_name];
    
    tic;
    NEV = openNEVCervical(fullfile(data_dir,[file_name '.nev']),'nosave');
    ns = openNSxCervical(fullfile(data_dir,[file_name '.ns5']));
    % extract spikes
    
    [ts,elec,wf] = deal([]);
    for i = 1:n_chans
        temp = findSpikesCervical(ns,'threshold',chan_thresholds(i),'channels',i,'filter',fc);
        ts = cat(2,ts,temp.TimeStamp');
        elec = cat(2,elec,temp.Electrode');
        wf = cat(2,wf, temp.Waveform);
    end
    [~,I] = sort(ts);
    spikes.TimeStamp = ts(I);
    spikes.Electrode = elec(I);
    spikes.Unit = [];
    spikes.Waveform = wf(:,I);
    
    disp('Found spikes. Saving...')
    
    NEV.Data.Spikes = spikes;
    
    % save spikes to file as NEV
    saveNEVSpikesCervical(NEV,data_dir,[file_name_save '.nev']);
    toc
end





