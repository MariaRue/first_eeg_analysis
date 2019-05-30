function [rt_fit,stats,rt_fit_all, stats_all] = rt_reward_history(data)




% empty vectors for regression for combined sessions
X_all = [];  
Y_all = [];  

num_sessions = length(data);

for i = 1:num_sessions
    
    % get rts 
    
    % for turtle eye - only use trials up to 1000 
    
    if length(data(i).data(:,1)) > 1000 
       data(i).data = data(i).data(1:1000,:);
     
    end 
%  
%     
    % normalise rts  = check distributions to make sure correct normal.
    % function is used either log or 1/rt should work 
  
    rt_norm = log(data(i).data(:,3)); %  1/rt or log(rts)
    
    sd = std(rt_norm);
    
    
    idx_rt_keep = rt_norm < (2.5 * sd);
    
    data_keep(i).data = data(i).data(idx_rt_keep,:);
    
    
%     a = num2str(i);
%     
%     figure
%     hist(data(i).data(:,3));
%     title(['untransformed',a])
%     
%     
%     
%     figure
%     hist(rt_norm); 
%     title(['transformed session',a])

% figure 
% hist(data(i).data(idx_rt_keep,3));
% title('rt keep')
    

% glm


% 
num_trials(i) = length(data_keep(i).data(:,1));  


% regressors for within session variance - subdivide each session into 2
% parts 
% num_trials05 = ceil(num_trials(i)./2);
% 
% session_part_1 =  [ones(num_trials05,1);zeros(num_trials(i)-num_trials05,1)];
% session_part_2 = [zeros(num_trials05,1); ones(num_trials(i)-num_trials05,1)];


% coherence = (data_keep(i).data(:,4) - mean(data_keep(i).data(:,4)))./std(data_keep(i).data(:,4));  

winright =  [0; (data_keep(i).data(1:num_trials(i)-1,5)== 1 & data_keep(i).data(1:num_trials(i)-1,2) == 1)]; % previous trial won right
winleft = [0; (data_keep(i).data(1:end-1,5)== 1 & data_keep(i).data(1:end-1,2) == 0)] .* -1; % previous trial won left%   

won_staircase2 = [0; (data_keep(i).data(1:num_trials(i)-1,5)== 1 & data_keep(i).data(2:num_trials(i),5) == 1)]; % last 2 trials won
won_staircase3 = [0; 0; (data_keep(i).data(1:num_trials(i)-2,5)== 1 & data_keep(i).data(3:num_trials(i),5) == 1)];  % last 3 trials won  

X{i} = [ones(num_trials(i),1),data_keep(i).data(:,4), winright, winleft]; % regressors
Y{i} = rt_norm(idx_rt_keep); % data

[rt_fit{i},~,stats{i}] = glmfit(X{i},Y{i},'normal','link','identity','constant','off'); %normal fit parameters


% save regressors and data for glm for combined session rts 
X_all = [X_all; X{i}]; 
Y_all = [Y_all; Y{i}];





end %loop through sessions 


% glm for combined sessions 

const_bet_sessions = [];

const_bet_sessions = [const_bet_sessions,[ones(num_trials(1),1); zeros(sum(num_trials(2:end)),1)]];



for i = 2:num_sessions
    
  
 const_bet_sessions = [const_bet_sessions,[zeros(sum(num_trials(1:i-1)),1);ones(num_trials(i),1);zeros(sum(num_trials(i+1:end)),1)]];
   
    
end

% X_all = [const_bet_sessions, X_all(:,2:end)];
% X_all = [const_bet_sessions, X_all];

[rt_fit_all,~,stats_all] = glmfit(X_all,Y_all,'normal','link','identity','constant','off'); %normal fit parameters
% 



end % end function