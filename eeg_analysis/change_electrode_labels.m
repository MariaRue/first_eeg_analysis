function [easy_cap_labels] = change_electrode_labels(labels)
% this function changes the electrod labels of the eeglab/curry files to
% to the same names used for the EEG cap layout from easy cap (EasycapM1.lay) we are using
% in fieldtrip to plot topo plots - rename from FP1 to Fp1 and CZ to Cz

% loop through all labels on cap (61 in our case) and change 2nd letter to
% lower case

easy_cap_labels = labels; 
for i  = 1:61 
    
    if i == 40 || i == 22 || i == 56 || i == 6
    
       easy_cap_labels{i}(3) = lower(labels{i}(3));   
   
    elseif i == 1 ||  i == 3 || i == 31 || i == 49 || i == 60 || i == 13
         easy_cap_labels{i}(2) = lower(labels{i}(2)); 
         
    elseif i == 2
        
        easy_cap_labels{i}(3) = lower(labels{i}(3));   
         easy_cap_labels{i}(2) = lower(labels{i}(2)); 
    end
    
    
end 



end 