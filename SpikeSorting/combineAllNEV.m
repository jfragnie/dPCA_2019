clear all;
close all;
clc;

spike_chans = 1:128;

%% Set up the location of the data
%   Need to define where the data is kept on this computer, which monkey
% you want, and which date you want
% root_dir = '/Users/mattperich/Dropbox/Research/Data/';
root_dir = 'C:\Users\jules\Downloads\Projet_Bach\HD_data\DataHD20190212';
monkey = 'Rey';
the_date = '20190220';
file_prefix = ['ReTh_' the_date '_' monkey '_Brain'];

% monkey = 'Jo';
% the_date = '20180615';
% file_prefix = ['ReTh_' the_date '_' monkey '_RecordBrain'];


%% set up the necessary info
data_dir = [fullfile(root_dir, monkey, the_date, 'BLACKROCK') filesep];

%% find all of the available .ns5 files
all_files  = dir(fullfile(data_dir,'*.ns5'));



%% loop along all files and rethreshold
all_data = [];
for iFile = 2:length(all_files)
    disp(['File ' num2str(iFile) ' of ' num2str(length(all_files))]);
    file_name = all_files(iFile).name(1:end-4);
    file_name_save = ['ReTh_' file_name];
    
tic
    ns = openNSxCervical(fullfile(data_dir,[file_name '.ns5']));

    all_data = [all_data, ns.Data(spike_chans,:)];
    toc
end

tic
writemda(all_data,fullfile(root_dir,monkey,the_date,'BLACKROCK',[monkey '_' the_date '_MountainSortData']),'int16');
toc


