function [logit_stats, logit_statsall] = PMF_fit_prev_trials(data,task)


% PMF fit with t-1 to t-4 trials won left/right 

num_sessions = length(data); %number of sessions
dataall = [];

for i = 1:num_sessions % loop through sessions
    
    historyr = zeros(length(data(i).data(:,1)),4);
    historyl = zeros(length(data(i).data(:,1)),4); 
    dataall = [dataall; data(i).data];
    
    
%     if task == 'c' 
%         
%         data = transform_touch_data(data);
%         
%         keep = data(i).data(:,4) ~= 0 ; 
%         
%         data(i).data = data(i).data(keep,:);
%         
%         
%     end
    
    num_trials = length(data(i).data(:,1)); % number of trials in session
    
%     winrightt1 = [(data(i).data(2:num_trials-1,5)== 1 & data(i).data(2:num_trials-1,2) == 1)]; % previous trial won right
%     winleftt1 = [(data(i).data(2:end-1,5)== 1 & data(i).data(2:end-1,2) == 0)]; % previous trial won left
    
    historyr(:,1) =  [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(1:num_trials-1,2) == 1)]; % previous trial won right
    historyl(:,1)  = [0; (data(i).data(1:end-1,5)== 1 & data(i).data(1:end-1,2) == 0)].*-1; % previous trial won left%   

    historyr(:,2) = [0; 0; (data(i).data(1:num_trials-2,5)== 1 & data(i).data(1:num_trials-2,2) == 1)]; % t-2 won right
    historyl(:,2) = [0; 0; (data(i).data(1:end-2,5)== 1 & data(i).data(1:end-2,2) == 0)].*-1; % t-2 won left
    
    historyr(:,3) = [0; 0; 0; (data(i).data(1:num_trials-3,5)== 1 & data(i).data(1:num_trials-3,2) == 1)]; % t-3 won right
    historyl(:,3) = [0; 0; 0; (data(i).data(1:end-3,5)== 1 & data(i).data(1:end-3,2) == 0)].*-1; % t-3 won left
    
    historyr(:,4) = [0; 0; 0; 0; (data(i).data(1:num_trials-4,5)== 1 & data(i).data(1:num_trials-4,2) == 1)]; % t-4 won right
    historyl(:,4) = [0; 0; 0; 0; (data(i).data(1:end-4,5)== 1 & data(i).data(1:end-4,2) == 0)].*-1; % t-4 won left
    
%     winright = [0; 0;  winrightt2 == 1];
%     winleft = [0; 0; winleftt2 == 1];
    
%     new_data_right = [new_data_right; data(i).data(logical(winright),:)]; % save trials corresponding to choice history to new matrix for plottting later on
%     new_data = [new_data_left; data(i).data(logical(winleft),:)];
%     
    Y{i} = [data(i).data(:,2)]; % column indicating choice
    X{i} = [ones(length(data(i).data(:,2)),1), data(i).data(:,4),historyr,historyl];
    
%     
%     % plot single sessions
%     
    [PMF_fit{i},~,logit_stats{i}] = glmfit(X{i},Y{i},'binomial','link','logit','constant','off'); %logistic fit parameters
    
%  
%     
%     
%     % get actual data points
%     [percent_right,~,confidence_interval,~,coherences] = PMF_data(data(i).data);
%     
%     % define x range over which fitted PMF is calculated
%     if task == 'd'
%         xrange = min(min(data(i).data(:,4)))-0.01:0.01:max(max(data(i).data(:,4)))+0.01;
%     elseif task == 'c'
%         xrange = min(min(data(i).data(:,4)))-10:0.01:max(max(data(i).data(:,4)))+10;
%     end
%     
%     % calculate fitted PMF curve
%     
%     param.b0 = PMF_fit{i}(1);
%     param.b1 = PMF_fit{i}(2);
%     
% 
%     
%     [p] = logist_PMF(xrange,param,0);
%  
%     
%     
%     
%     % plot fitted PMF with actual data
%     figure
%     hold on
%     plot(xrange,p.*100,'k-'); % fitted PMF
% 
%     errorbar(coherences,percent_right,confidence_interval,'ko') % data
%     
%     hold off
%     ylim([0 100]);
%     
%     if task == 'd'
%         xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
%     elseif task == 'c'
%         xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
%     end
%     title(sprintf( 'PMFs choice history session %d', i), 'FontSize', 14);
%     ylabel('percentage choice to the right', 'FontSize', 14)
    
end
% combined data

X_all = vertcat(X{:});
Y_all = vertcat(Y{:});



[PMF_fitall,~,logit_statsall] = glmfit(X_all,Y_all,'binomial','link','logit','constant','off'); %logistic fit parameters




% get actual data points
[percent_right_all,~,confidence_interval_all,~,coherences] = PMF_data(dataall);

% 
% % define x range over which fitted PMF is calculated
% if task == 'd'
%     xrange = min(min(new_data_left(:,4)))-0.01:0.01:max(max(new_data_left(:,4)))+0.01;
% elseif task == 'c'
%     xrange = min(min(new_data_left(:,4)))-10:0.01:max(max(new_data_left(:,4)))+10;
% end

% % calculate fitted PMF curve
% 
% param.b0 = PMF_fitall(1);
% param.b1 = PMF_fitall(2);
% 
% 
% 
% [p_all] = logist_PMF(xrange,param,0);
% 
% 
% 
% 
% % plot fitted PMF with actual data
% figure
% hold on
% plot(xrange,p_all.*100,'k-'); % fitted PMF
% errorbar(coherences,percent_right_all,confidence_interval_all,'ko') % datahold off
% ylim([0 100]);
% 
% 
% if task == 'd'
%     xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
% elseif task == 'c'
%     xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
% end
% title( 'combined data', 'FontSize', 14);
% ylabel('percentage choice to the right', 'FontSize', 14)
% 
% 


end