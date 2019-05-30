function [logit_stats] = PMF_history_even_num(data,task)

% calculates logistic PMF fit taking into account reward history of
% previous trial in form of two vectors - the winning vector with correct
% right = 1, correct left, -1, and false 0 and a loosing vector with false
% right 1, false left -1 and correct = 0 
%
% Input is behavioural data from TTT_training_analysis. For touch screen
% data, data has to be transformed with transform_touch_data.m. 

% Outputs are PMF fits either for each session or for combined behavioural
% data 

% Maria Ruesseler, University of Oxford, 2017 


if isstruct(data) % session data 
    
    num_sessions = length(data); % number of sessions 
    PMF_fit = zeros(num_sessions,4); % set up vector to save PMF fit parmaters
    
    for j = 1 : num_sessions % loop through sessions 
        
        num_trials = size(data(j).data(:,1),1); % get number of trials in session
        winning = zeros(floor(num_trials/2),1); % set up winning vector 
        loosing = zeros(floor(num_trials/2),1); % set up loosign vector 
        
        
        
        z = 1;
        for i = 2 :2 : num_trials % loop through trials 
            
            if data(j).data(i-1,5) == 1 && data(j).data(i-1,2) == 1 % correct rightward choice
                
                winning(z) = 1;
                
            elseif data(j).data(i-1,5) == 1 && data(j).data(i-1,2) == 0 % correct leftward choice
                
                winning(z) = -1;
                
            elseif data(j).data(i-1,5) == 0 % false choice
                
                winning(z) = 0;
                
            end % if winning
            
            
            
            if data(j).data(i-1,5) == 0 && data(j).data(i-1,2) == 1 % failed rightward choice
                
                loosing(z) = 1;
                
            elseif data(j).data(i-1,5) == 0 && data(j).data(i-1,2) == 0 % failed leftward choice
                
                loosing(z) = -1;
                
            elseif data(j).data(i-1,5) == 1 % correct choice
                
                loosing(z) = 0;
                
            end % if loosing
            
            z = z + 1;
            
        end % loop through trials in session
        
        

        winleft = [winning(:,:).*(1-data(j).data(1:2:end-1,2))]; %won on left on t-1?
        winright = [winning(:,:).*data(j).data(1:2:end-1,2)]; % won on right on t-1?
        choseright = [0; data(j).data(1:2:end-1,2)]; %chose right on t-1?
        
        Y = data(j).data(2:2:end,2); % column indicating rightward choice
        
        X = [ones((floor(num_trials/2)),1),data(j).data(2:2:end,4), winleft, winright];
        %first column is constant, second column list with coherences, 3rd column
        %winning vector, 4th loosing vector
        
        
        [PMF_fit(j,:),~,logit_stats{j}] = glmfit(X,Y,'binomial','link','logit','constant','off'); %logistic fit parameters
        
        
        param.b0 =  PMF_fit(j,1); % constant
        param.b1 =  PMF_fit(j,2); % coherence
        param.b2 =  PMF_fit(j,3); % winning 
        param.b3 =  PMF_fit(j,4); % loosing 
        
        
        
        
        
        % get actual data points 
        [percent_right,n,confidence_interval,rt,coherences] = PMF_data(data(j).data(2:2:end,:));
        
        
        % define x range over which fitted PMF is calculated 
        if task == 'd'
            xrange = min(min(data(j).data(:,4)))-0.01:0.01:max(max(data(j).data(:,4)))+0.01;
        elseif task == 'c'
            xrange = min(min(data(j).data(:,4)))-10:0.01:max(max(data(j).data(:,4)))+10;
        end
        
        % calculate fitted PMF curve 
        [p] = logist_PMF(xrange,param,0);
        
        % plot fitted PMF with actual data 
        subplot(ceil(num_sessions/2),2,j)
        hold on
        plot(xrange,p.*100,'g-'); % fitted PMF
        errorbar(coherences,percent_right,confidence_interval,'go') % data
        hold off
        ylim([0 100]);
        
        if task == 'd'
            xlabel('CW rotation          disparity        CCW rotation');
        elseif task == 'c'
            xlabel('leftward motion       coherence     rightward motion');
        end
        title(sprintf('Data for session %d. b0 = %.2f; b1 = %.2f. b2 = %.2f b3 = %.2f', j, PMF_fit(j,1), PMF_fit(j,2), PMF_fit(j,3), PMF_fit(j,4)), 'FontSize', 14);
        
        
        
        
        
        
    end % loop through sessions
    
    
else % combined data
    
    
    num_trials = size(data(:,1),1);
    winning = zeros(num_trials,1);
    loosing = zeros(num_trials,1);
    
    
    
    for i = 2 : num_trials
        
        if data(i-1,5) == 1 && data(i-1,2) == 1 % correct rightward choice
            
            winning(i) = 1;
            
        elseif data(i-1,5) == 1 && data(i-1,2) == 0 % correct leftward choice
            
            winning(i) = -1;
            
        elseif data(i-1,5) == 0 % false choice
            
            winning(i) = 0;
            
        end
        
        if data(i-1,5) == 0 && data(i-1,2) == 1 % failed rightward choice
            
            loosing(i) = 1;
            
        elseif data(i-1,5) == 0 && data(i-1,2) == 0 % failed leftward choice
            
            loosing(i) = -1;
            
        elseif data(i-1,5) == 1 % correct choice
            
            loosing(i) = 0;
            
        end
        
    end
    
    
    Y = data(:,2); % column indicating rightward choice
    
    X = [ones(num_trials,1),data(:,4),winning,loosing];
    %first column is constant, second column list with coherences, 3rd column
    %indicates whether previous trial correct and whether it was a left (-1)
    %or rightward (+1) choice, otherwise false choice = 0
    
    [PMF_fit,~,logit_stats] = glmfit(X,Y,'binomial','link','logit','constant','off'); %logistic fit parameters
    
    
    param.b0 =  PMF_fit(1);
    param.b1 =  PMF_fit(2);
    param.b2 =  PMF_fit(3);
    param.b3 =  PMF_fit(4);
    
   
    
    [percent_right,n,confidence_interval,rt,coherences] = PMF_data(data);
    
    
    if task == 'd'
        xrange = min(min(data(:,4)))-0.01:0.01:max(max(data(:,4)))+0.01;
    elseif task == 'c'
        xrange = min(min(data(:,4)))-10:0.01:max(max(data(:,4)))+10;
    end
    
    
    
    [p] = logist_PMF(xrange,param,0);
    
    figure
    hold on
    plot(xrange,p.*100,'k-');
    errorbar(coherences,percent_right,confidence_interval,'ko')
    hold off
    if task == 'd'
        xlabel('CW rotation          disparity        CCW rotation');
    elseif task == 'c'
        xlabel('leftward motion       coherence     rightward motion');
    end
    
    
    
    title(sprintf('combined data. b0 = %.2f; b1 = %.2f. b2 = %.2f b3 = %.2f', PMF_fit(1), PMF_fit(2), PMF_fit(3), PMF_fit(4)), 'FontSize', 14);
    
    ylim([0 100]);
    

end




end