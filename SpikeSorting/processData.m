function [data_mod,Type] = processData(data)
%Sorting the data by electrode type. FinalData.ElectrodeLabel is either
%elec1-#, elec2-#, elec3-#

% FIRST IDEA
%creating a matrix of 1 and zero which defines with elec type.
ElecType1_v = cell(1,128);
for i=1:length(data)
   ElecType1 = data(i).ElectrodeLabel;
   ElecType1 = ElecType1';  %convertCharsToStrings(A');
   ElecType1_v{i} = ElecType1;
end

%Vec type. The first column represent all the vector of type elec1 (if val=1)
Type = zeros(128,3);
electrodes= {'elec1','elec2','elec3'};
for i=1:3
     Rep = ~cellfun('isempty',strfind(ElecType1_v,electrodes{i}));
     Type(:,i)= double(Rep');
end

%SECOND IDEA
% Already change the order of the data s.t. elec1 are bottom followed by
% elec2 and then elec3 such that data is ordered by electrode type.
elec1 = struct([]); elec2 = struct([]); elec3 = struct([]);

for i=1:length(data)
   ElecVal= data(i).ElectrodeLabel;
   ElecVal = cellstr(ElecVal');
   if(~cellfun('isempty',strfind(ElecVal,electrodes{1})))
        elec1 = [data(i) elec1];   
   elseif(~cellfun('isempty',strfind(ElecVal,electrodes{2})))
        elec2 = [data(i) elec2];
   else
        elec3 = [data(i) elec3];
   end
end
data_mod = [elec1 elec2 elec3];     
clear elec1 elec2 elec3 i Rep ElecVal ElecType1 ElecType1_v electrodes;
end

