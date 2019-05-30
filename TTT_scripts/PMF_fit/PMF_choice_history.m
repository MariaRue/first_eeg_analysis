function [PMF_fit,logit_stats,PMF_fit_all, logit_statsall] = ...
    PMF_choice_history(data,task)


num_sessions = length(data); %number of sessions

for i = 1:num_sessions % loop through sessions
    
    if task == 'c'
        
        data = transform_touch_data(data);
        
        keep = data(i).data(:,4) ~= 0;
        
        data(i).data = data(i).data(keep,:);
        
        
    end
    
    %keep = data(i).data(:,4) >  - 0.02 & data(i).data(:,4) < 0.02;
    
    %data(i).data = data(i).data(keep,:);
    
    num_trials = length(data(i).data(:,1)); % number of trials in session
    
    winright{i} = [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(1:num_trials-1,2) == 1)]; % previous trial won right
    winleft{i}  = [0; (data(i).data(1:end-1,5)== 1 & data(i).data(1:end-1,2) == 0)].*-1; % previous trial won left
    
    
    
    %
    %     winrightt1 = [(data(i).data(2:num_trials-1,5)== 1 & data(i).data(2:num_trials-1,2) == 1)]; % previous trial won right
    %     winleftt1 = [(data(i).data(2:end-1,5)== 1 & data(i).data(2:end-1,2) == 0)]; % previous trial won left
    %
    %     winrightt2 = [(data(i).data(1:num_trials-2,5)== 1 & data(i).data(1:num_trials-2,2) == 1)]; % t-2 won right
    %     winleftt2 = [(data(i).data(1:end-2,5)== 1 & data(i).data(1:end-2,2) == 0)]; % t-2 won left
    %
    %     winright{i} = [0; 0; winrightt1 == 1 & winrightt2 == 1];
    %     winleft{i} = [0; 0; winleftt1 == 1 & winleftt2 == 1].* -1;
    
    Y{i} = data(i).data(:,2); % column indicating rightward choice
    X{i} = [ones(num_trials,1),data(i).data(:,4),winright{i}, winleft{i}];
    
    
    
    % plot single sessions
    
    [PMF_fit{i},~,logit_stats{i}] = glmfit(X{i},Y{i},'binomial','link','logit','constant','off'); %logistic fit parameters
    
    % get actual data points
    [percent_right,~,confidence_interval,~,coherences] = PMF_data(data(i).data);
    
    
    % define x range over which fitted PMF is calculated
    if task == 'd'
        xrange = min(min(data(i).data(:,4)))-0.01:0.01:max(max(data(i).data(:,4)))+0.01;
    elseif task == 'c'
        xrange = min(min(data(i).data(:,4)))-10:0.01:max(max(data(i).data(:,4)))+10;
    end
    
    % calculate fitted PMF curve
    
    param.b0 = PMF_fit{i}(1);
    param.b1 = PMF_fit{i}(2);
    
    
    [p] = logist_PMF(xrange,param,0);
    
    
    
    
    % plot fitted PMF with actual data
    figure
    hold on
    plot(xrange,p.*100,'k-'); % fitted PMF
    
    errorbar(coherences,percent_right,confidence_interval,'ko') % data
    
    hold off
    ylim([0 100]);
    
    
    if task == 'd'
        xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
    elseif task == 'c'
        xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
    end
    title(sprintf( 'PMFs choice history session %d', i), 'FontSize', 14);
    ylabel('percentage choice to the right', 'FontSize', 14)
    
    
    
end



% combined data

Xall = vertcat(X{:});

Yall = vertcat(Y{:});


[PMF_fit_all,~,logit_statsall] = glmfit(Xall,Yall,'binomial','link','logit','constant','off'); %logistic fit parameters


all_data = vertcat(data(:).data);

% get actual data points
[percent_rightall,~,confidence_intervalall,~,coherences] = PMF_data(all_data);


% define x range over which fitted PMF is calculated
if task == 'd'
    xrange = min(min(all_data(:,4)))-0.01:0.01:max(max(all_data(:,4)))+0.01;
elseif task == 'c'
    xrange = min(min(all_data(:,4)))-10:0.01:max(max(all_data(:,4)))+10;
end

% calculate fitted PMF curve

paramall.b0 = PMF_fit_all(1);
paramall.b1 = PMF_fit_all(2);



[pall] = logist_PMF(xrange,paramall,0);



% plot fitted PMF with actual data
figure
hold on
plot(xrange,pall.*100,'k-'); % fitted PMF

errorbar(coherences,percent_rightall,confidence_intervalall,'ko') % data

hold off
ylim([0 100]);



if task == 'd'
    xlabel('CW rotation          disparity        CCW rotation', 'FontSize', 14);
elseif task == 'c'
    xlabel('leftward motion       coherence     rightward motion', 'FontSize', 14);
end
title( 'combined data', 'FontSize', 14);
ylabel('percentage choice to the right', 'FontSize', 14)








end