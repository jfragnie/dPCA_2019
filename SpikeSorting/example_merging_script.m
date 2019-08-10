clear all;
close all;
clc;

%% Set up the location of the data
%   Need to define where the data is kept on this computer, which monkey
% you want, and which date you want
% root_dir = '/Users/mattperich/Dropbox/Research/Data/';
%root_dir = 'C:\Users\Matt\Dropbox\Research\Data\';
root_dir = 'C:\Users\jules\Downloads\Projet_Bach\HD_data\';

% monkey = 'Cersei';
% the_date = '20180615';
% file_prefix = ['ReTh_' the_date '_' monkey '_Brain'];

monkey = 'Rey';
the_date = '20190212';
file_prefix = ['ReTh_' the_date monkey];% '_' monkey '_RecordBrain'];


%% set up the necessary info
data_dir = [fullfile(root_dir, monkey, the_date, 'BLACKROCK_WIRED') filesep];


%% First, rethreshold all the .NS5 files and save them as new .NEVs
%   Each .NS5 file in the destination folder will be rethresholded and
% coverted to a new .NEV file with 'ReTh_' added to the beginning.
%   e.g., 'ReTh_Sansa_20180924_Brain001.nev'
rethresholdAllNEV;


%% Second, merge all of the new .NEVs into a single .NEV
%   This will save a merged version of all of the ReTh_ nev files with
% '-spikes' appended to the end of the root filename.
%   e.g., 'ReTh_Sansa_20180924_Brain-spikes.nev'
tic;
mergingStatus = mattProcessSpikesForSorting(data_dir,file_prefix,false);
toc;

%% Third, sort all of the spikes
%   Do this in offline sorter. To make the next part work, save the file
% with the same filename but add '-s'
%   e.g., 'Sansa_20180924_Brain-spikes-s.nev'


%% Fourth, separate the sorted file back into their own files
%   This will split the sorted .NEV into a .mat file representing each of
% the original ReTh_ NEV files. The .mat file will have the same name but
% contain the NEV struct that would be loaded by a function call such as
% openNEV.

%sorted them first!
mergingStatus = mattProcessSpikesForSorting(data_dir,file_prefix,true);


