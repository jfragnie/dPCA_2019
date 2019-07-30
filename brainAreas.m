function [AVG1,AVG2,AVG3] = brainAreas(Type,averageFR)
%Separation of data averaged in three diff brain areas

AVG1 = [];  AVG2 = [];  AVG3 = [];
for i=1:length(Type)
   if (Type(i,1) == 1)
        AVG1 = [averageFR(i,:); AVG1] ;
   elseif (Type(i,2) == 1)
        AVG2 = [averageFR(i,:); AVG2]; 
   else 
        AVG3 = [averageFR(i,:); AVG3]; 
   end
end

end

