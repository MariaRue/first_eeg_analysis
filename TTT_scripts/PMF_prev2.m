function [logit_stats_right, logit_stats_left, logit_stats_rightall, logit_stats_leftall] = PMF_prev2(data,task)

% calculates same as PMF_prev1 just for t-2. 


num_sessions = length(data); %number of sessions
new_data_right = [];
new_data_left = [];

for i = 1:num_sessions % loop through sessions
    
    
%     
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
    
    winrightt1 = [(data(i).data(2:num_trials-1,5)== 1 & data(i).data(2:num_trials-1,2) == 1)]; % previous trial won right
    winleftt1 = [(data(i).data(2:end-1,5)== 1 & data(i).data(2:end-1,2) == 0)]; % previous trial won left
    
    winrightt2 = [(data(i).data(1:num_trials-2,5)== 1 & data(i).data(1:num_trials-2,2) == 1)]; % t-2 won right
    winleftt2 = [(data(i).data(1:end-2,5)== 1 & data(i).data(1:end-2,2) == 0)]; % t-2 won left
    
    winright = [0; 0;  winrightt2 == 1];
    winleft = [0; 0; winleftt2 == 1];
    
    new_data_right = [new_data_right; data(i).data(logical(winright),:)]; % save trials corresponding to choice history to new matrix for plottting later on
    new_data_left  = [new_data_left; data(i).data(logical(winleft),:)];
    
    Y_right{i} = [data(i).data(logical(winright),2)]; % column indicating rightward choice
    X_right{i} = [ones(sum(winright),1),data(i).data(logical(winright),4)];
    
    Y_left{i} = data(i).data(logical(winleft),2); % column indicating rightward choice
    X_left{i} = [ones(sum(winleft),1),data(i).data(logical(winleft),4)];
    
    % plot single sessions
    
    [PMF_fitright{i},~,logit_stats_right{i}] = glmfit(X_right{i},Y_right{i},'binomial','link','logit','constant','off'); %logistic fit parameters
    
    [PMF_fitleft{i},~,logit_stats_left{i}] = glmfit(X_left{i},Y_left{i},'binomial','link','logit','constant','off'); %logistic fit parameters
    
    
    % get actual data points
    [percent_right_r,~,confidence_interval_right,~,~] = PMF_data(data(i).data(logical(winright),:));
    [percent_right_l,~,confidence_interval_left,~,coherences] = PMF_data(data(i).data(logical(winleft),:));
    
    % define x range over which fitted PMF is calculated
    if task == 'd'
        xrange = min(min(new_data_left(:,4)))-0.01:0.01:max(max(new_data_left(:,4)))+0.01;
    elseif task == 'c'
        xrange = min(min(new_data_left(:,4)))-10:0.01:max(max(new_data_left(:,4)))+10;
    end
    
    % calculate fitted PMF curve
    
    paramright.b0 = PMF_fitright{i}(1);
    paramright.b1 = PMF_fitright{i}(2);
    
    paramleft.b0 = PMF_fitleft{i}(1);
    paramleft.b1 = PMF_fitleft{i}(2);
    
    [p_right] = logist_PMF(xrange,paramright,0);
    [p_left] = logist_PMF(xrange,paramleft,0);
    
    
    
    % plot fitted PMF with actual data
    figure
    hold on
    plot(xrange,p_right.*100,'k-'); % fitted PMF
    plot(xrange,p_left.*100,'b-'); % fitted PMF
    errorbar(coherences,percent_right_r,confidence_interval_right,'ko') % data
    errorbar(coherences,percent_right_l,confidence_interval_left,'bo') % data
    hold off
    ylim([0 100]);
    
    legend('t-1 won right', 't-1 won left')
    
    if task == 'd'
        xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
    elseif task == 'c'
        xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
    end
    title(sprintf( 'PMFs choice history session %d', i), 'FontSize', 14);
    ylabel('percentage choice to the right', 'FontSize', 14)
    
end
% combined data

X_rightall = vertcat(X_right{:});
Y_rightall = vertcat(Y_right{:});
X_leftall = vertcat(X_left{:});
Y_leftall = vertcat(Y_left{:});


[PMF_fitrightall,~,logit_stats_rightall] = glmfit(X_rightall,Y_rightall,'binomial','link','logit','constant','off'); %logistic fit parameters

[PMF_fitleftall,~,logit_stats_leftall] = glmfit(X_leftall,Y_leftall,'binomial','link','logit','constant','off'); %logistic fit parameters


% get actual data points
[percent_right_rall,~,confidence_interval_rightall,~,~] = PMF_data(new_data_right);
[percent_right_lall,~,confidence_interval_leftall,~,coherences] = PMF_data(new_data_left);

% define x range over which fitted PMF is calculated
if task == 'd'
    xrange = min(min(new_data_left(:,4)))-0.01:0.01:max(max(new_data_left(:,4)))+0.01;
elseif task == 'c'
    xrange = min(min(new_data_left(:,4)))-10:0.01:max(max(new_data_left(:,4)))+10;
end

% calculate fitted PMF curve

paramright.b0 = PMF_fitrightall(1);
paramright.b1 = PMF_fitrightall(2);

paramleft.b0 = PMF_fitleftall(1);
paramleft.b1 = PMF_fitleftall(2);

[p_rightall] = logist_PMF(xrange,paramright,0);
[p_leftall] = logist_PMF(xrange,paramleft,0);



% plot fitted PMF with actual data
figure
hold on
plot(xrange,p_rightall.*100,'k-'); % fitted PMF
plot(xrange,p_leftall.*100,'b-'); % fitted PMF
errorbar(coherences,percent_right_rall,confidence_interval_rightall,'ko') % data
errorbar(coherences,percent_right_lall,confidence_interval_leftall,'bo') % data
hold off
ylim([0 100]);

legend('t-1 won right', 't-1 won left')

if task == 'd'
    xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
elseif task == 'c'
    xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
end
title( 'combined data', 'FontSize', 14);
ylabel('percentage choice to the right', 'FontSize', 14)


end