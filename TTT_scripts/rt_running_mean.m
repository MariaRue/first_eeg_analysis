%% script 
% requires data struct with single sessions and calculates the running mean
% of rt across trials for each session 

num_sessions  = length(data); 
figure;
for i = 1 :num_sessions
    
    
    rts = data(i).data(:,3);
    num_trials = length(rts);
    
    
    running_mean{i} = zeros((num_trials-20));
    j = 1; % trial counter
    
    while j <= (num_trials-20)
        
        running_mean{i}(j) = mean(rts(j:j+19));
       
        j = j+1;
    end %loop through trials 
  
    hold on 
    plot([1:num_trials-20],running_mean{i})
    
    
end % loop through sessions 

xlabel('trials')
ylabel('running mean (s)')
title('running mean')
hold off
