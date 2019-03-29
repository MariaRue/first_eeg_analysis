function [upsampled_triggers] = upsample_triggers(triggers_to_upsample,fs_old, fs_new,variable)  

% function to upsample behav trigger vals to match behav triggers and
% eyetracker triggers 
% Input: triggers_to_upsample = vector with triggers 
%        fs_old               = hz in which triggers_to_upsample was
%                               sampled
%        fs_new               = hz to which triggers should be upsampled to
%        



    
  
nsamples_insert = round(fs_new/fs_old); % num of samples between old and new trigger postion
upsampled_triggers = zeros(length(triggers_to_upsample).*nsamples_insert,1); % size of vector with upsampled triggers


switch variable 
    case 'trigger'

upsampled_triggers(1:nsamples_insert:end) = triggers_to_upsample; % triggers at new fs 
    case 'coherence'
for i = 1:length(triggers_to_upsample)
if i ==1 
    
    start_trig = 1;
    
end 

end_trig = start_trig + nsamples_insert-1;

    upsampled_triggers(start_trig:end_trig) = ones(length(start_trig:end_trig),1) .* triggers_to_upsample(i); 
    
    start_trig = end_trig + 1; 
    
end

end





end 