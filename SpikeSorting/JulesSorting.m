%script to prepare the matrix:
%the matrix is a units x timePoint matrix.
%If there is a unit. the matrix will contain 1, either filled with zeros.

root_dir = 'C:\Users\jules\Downloads\Projet_Bach\HD_data\';
monkey = 'Rey';

%% 12.02.19

the_date = '20190212';
file_path = [fullfile(root_dir, monkey, the_date, 'BLACKROCK_WIRED') filesep];
file_path2 = [fullfile(root_dir, monkey, the_date, 'VICON') filesep]; 
NEV_20190212 = [];
data_mod = {};
unitsPerElec1 = zeros(128,24);
%unitsTime_v = {};
%whichElec_v = {};
%combOverTrials_v = {};

for i=1
    if (i < 10)
        file_prefix = ['ReTh_', the_date, monkey '-spikes-s'];
        file_prefix = ['ReTh_', the_date, monkey 'NSP100' num2str(i)];
        file_prefix2 = [the_date, '_', monkey, '_Training0' num2str(i)];
        file_prefix3 = [the_date, '_', monkey, '_Training0' num2str(i) '_indMov'];
    else 
        file_prefix = ['ReTh_' the_date monkey 'NSP10' num2str(i)];
        file_prefix2 = [the_date, '_', monkey, '_Training' num2str(i)];
        file_prefix3 = [the_date, '_', monkey, '_Training' num2str(i) '_indMov'];
    end
load([file_path file_prefix]);
load([file_path2 file_prefix2]);
load([file_path2 file_prefix3]);
    for elec=1:128
        prov = find(NEV.Data.Spikes.Electrode==elec);
        unit = max(NEV.Data.Spikes.Unit(prov));
        if ~isempty(unit)
            unitsPerElec1(elec,i) = unit;
        end
    end
    %Now sort the electrodes per brain areas (as I did in KOBAK folder)
    NEV.ElectrodesInfo = NEV.ElectrodesInfo(1:128);
    [data_mod,~] = processData(NEV.ElectrodesInfo);
    NEV.ElectrodesInfo = data_mod;
    %units = sum([NEV.Data.Spikes.Unit]); way too big as rep. of same elec
    units = sum(unitsPerElec1(:,i)); % equals 153 modulo 1 in this case..
    unitsTime = zeros(units,NEV.MetaTags.DataDuration); 
    
    %struct to know to which units correspond which electrodes!
    %the value in each 'case' is the corresponding electrodes!
    whichElec=zeros(units,1);
    old_nb = 0;
    for j=1:128
        new_nb = old_nb + double(NEV.ElectrodesInfo(j).Units);
        whichElec(old_nb+1:new_nb,1) = NEV.ElectrodesInfo(j).ElectrodeID;
        old_nb = new_nb;
    end
    
    DataLength = 1:NEV.MetaTags.DataDuration;
    for time=DataLength
        ind = find(NEV.Data.Spikes.TimeStamp==time);
        if ~isempty(ind)
            %units corresponding to the indices found in time stamps.
            unitt = NEV.Data.Spikes.Unit(ind);
            if ~isempty(unitt)
                %electrodes corresponding to the units found.
                elecUnit = NEV.Data.Spikes.Electrode(ind);
                
                %we have to compute the pos for each elec!
                %it will be the sum of all unit before it.
                for j=1:length(elecUnit)
                    %cannot sum it the same way is the electrode as
                    %multiple units
                    if (unitsPerElec1(elecUnit(j))==1)
                        pos(j) = sum(unitsPerElec1(1:elecUnit(j)));
                        
                    else
                        pos(j) = sum(unitsPerElec1(1:elecUnit(j)-1)) +1;
                        %the pos of the first el has to be @pos 1 one the array
                    end
                    unitsTime(pos,ind) = 1;
                end
                pos = [];
            end
        end
    end
        
    H = TimeStamp_finder(analog,kinematic,NEV,indMovAdjBeg);
    combOverTrials = TimeCutSort(H,unitsTime,NEV); 
    combOverTrials_v{i} = combOverTrials;    
    %unitsTime_v{i} = unitsTime;
    %whichElec_v{i} = whichElec;
    clear NEV unitsTime whichElec;
end

%% 20.02.19

% the_date = '20190220';
% file_path = [fullfile(root_dir, monkey, the_date, 'BLACKROCK_WIRED') filesep];
% NEV_20190220 = [];
% 
% for i=1:23
%     if (i < 10)
%         file_prefix = ['ReTh_' the_date monkey 'InjuryNSP100' num2str(i)];
%     else 
%         file_prefix = ['ReTh_' the_date monkey 'InjuryNSP10' num2str(i)];
%     end
% load([file_path file_prefix]);
% NEV_20190220 = [NEV NEV_20190220];
% clear NEV ;
% end


%% clearer
%clear elec elecUnit event file_path file_prefix H 
%clear i iA ind ind_NSP1 j neg NEV_20190212 
%clear new_nb old_nb pos prov unit units unitt uuu



%% find FR
%choose a value for elec!
elec=1; %first elec (on 128)
prov = find(NEV.Data.Spikes.Electrode==elec);
timeElec = NEV.Data.Spikes.TimeStamp(prov);%/300;  %in second
diff = double(diff(timeElec));
inv = 1./diff;