function [fit,stats,fit_all,stats_all] = reward_staircase_history(data)


% glm with regrossors for reward staircase

X_all = [];  
Y_all = []; 



num_sessions = length(data);

for i = 1:num_sessions
    
num_trials = length(data(i).data(:,1));    
  
    winright = [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(1:num_trials-1,2) == 1)]; % previous trial won right
    winleft  = [0; (data(i).data(1:end-1,5)== 1 & data(i).data(1:end-1,2) == 0)] .* -1; % previous trial won left

won_staircase2 = [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(2:num_trials,5) == 1)]; % last 2 trials won
won_staircase3 = [0; 0; (data(i).data(1:num_trials-2,5)== 1 & data(i).data(3:num_trials,5) == 1)];  % last 3 trials won  
   


X{i} = [ones(num_trials,1),data(i).data(:,4), winright, winleft, won_staircase2,  won_staircase3]; % regressors
Y{i} = data(i).data(:,2); % data
    
[fit{i},~,stats{i}] = glmfit(X{i},Y{i},'binomial','link','logit','constant','off'); %normal fit parameters

% save regressors and data for glm for combined session rts 
X_all = [X_all; X{i}]; 
Y_all = [Y_all; Y{i}];



end % loop through sessions 

[fit_all,~,stats_all] = glmfit(X_all,Y_all,'binomial','link','logit','constant','off'); %normal fit parameter




end % end function