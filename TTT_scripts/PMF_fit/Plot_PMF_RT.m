function [PMF_fit] = Plot_PMF_RT(data,logs,Gauss,task)

% returns parameter fits (mean, sd for cummulative gaussian PMF fitting, b0, b1 for logistic) for each
% session or for combined session data and plots respective PMF as well as
% rts for correct and failed trials

% Input:
% data = either structure with field for each session or matrix for
% combined session data created with TTT_analysis_eye or TTT_analysis
% log or Gauss = logicals, 1 for the PMF fitting algorithm one wants to
% use, 0 the other one
% task = 'd' for disparity/cylinder, 'c' for rdk/Shadlen

% Maria Ruesseler, University of Oxford, 2017




if isstruct(data) % session data
    
    session_num = length(data); % number of sessions 
    
    PMF_fit = zeros(length(session_num),2); % matrix in which we save the fitting parameters 
    
    
    
    for j = 1:session_num %loop through sessions
        
        % transform rts to exclude rts that a longer than 2.5 sds away from
        % mean 
        
        % transform data into a gaussian distribution (for Turtle touch I
        % would recomment 1/rt)
        
   rt_norm = log(data(j).data(:,3));
    
   % get the standard deviation 
   sd = std(rt_norm);
   
   % indices of trials that are below 2 1/2 times of sd - these are the
   % indices of the trials we wnat to keep
   idx_rt_keep = rt_norm < (2.5 * sd);
   
   % get the data with rts below 2.5 * sd 
   data_keep = data(j).data(idx_rt_keep,:);
   
         % get data info for plotting Rts and PMFs
        [percent_right,n,confidence_interval,rt,coherences] = PMF_data(data_keep); % (info needed for Y for glmfit logistic fit)
        
        if task == 'd' % if task is cylinder
            xrange = min(min(data(j).data(:,4)))-0.01:0.01:max(max(data(j).data(:,4)))+0.01; % range over which we evaluate fitted PMF
        elseif task == 'c' % if task is shadlen 
            xrange = min(min(data(j).data(:,4)))-10:0.01:max(max(data(j).data(:,4)))+10;
            
        end
        
       
        
        % fit PMF with GAussian or Logistic regression 
        if Gauss
            
            % how to calculate initial guesses?
            x0 = [1,1]; % initial guesses for mean(x0(1)) and sd (x0(2))
            
            options = optimset(); % set options to default value - for fminsearch function
            
            [PMF_fit(j,:)] = fminsearch(@cum_Gauss_loglike,x0,options,data(j).data); % find mean/sd
            
            param.b0 = PMF_fit(j,1); %mean
            param.b1 = PMF_fit(j,2); %sd
            
            [p] = cum_Gauss_PMF(xrange,param,0); % calculate cummulative Gaussian values at each xrange value for plotting later on
            
        elseif logs % fit logistic
            
            Y = [percent_right .* n ./100, n]; % first column number of rightward (CCW) choices for each stim, second column number of stim representations
            
            X = [ones(length(coherences),1),coherences]; % first column is constant, second column list with coherences
            
            PMF_fit(j,:) = glmfit(X,Y,'binomial','link','logit','constant','off'); %logistic fit parameters
            
            % parameter input in logist_PMF function which calculates
            % values for each point in xrange for plotting
            param.b0 = PMF_fit(j,1); 
            param.b1 = PMF_fit(j,2);
           
            %
            [p] = logist_PMF(xrange,param,0); %calculate fitted PMF for each point in xrange
            
        end %if Gauss
        
        
        %plot PMFs
%         figure
figure (1)
subplot(6,2,j)
        hold on
        plot(xrange,p.*100,'k-'); %plot fitted pmf
        errorbar(coherences,percent_right,confidence_interval,'ko') %plot behav. data
        hold off
        ylim([0 100]);
        
        
        %plot labels depending on rdk or cylinder task
        if task == 'c' %rdk
            xlabel('(leftward motion)              coherence            (rigthward motion)')
            ylabel('percentage rightward motion')
            
            
        elseif task == 'd' %cylinder
            xlabel('(CW rotation)                  disparity              (CCW rotation)')
            ylabel('percentage CCW rotation')
            
        end %plot label
        
        % title for GAuss or log
        if Gauss
            title(sprintf('Data for session %d. Mean = %.2f; SD = %.2f.', j, PMF_fit(j,1), PMF_fit(j,2)), 'FontSize', 14);
            
        else %log
            
            title(sprintf('Data for session %d. b0 = %.2f; b1 = %.2f.', j, PMF_fit(j,1), PMF_fit(j,2)), 'FontSize', 14);
            
        end %if Gauss title
        
        
        
        %Plot Rts
        
%         figure % correct rts
figure (2)
subplot(6,2,j)
        errorbar(coherences,rt(:,1),rt(:,2),rt(:,5),'ko-')
        hold on
        errorbar(coherences,rt(:,3),rt(:,4),rt(:,6),'rd-')
        hold off
        
        %plot labels depending on rdk or cylinder task
        if task == 'c'
            xlabel('(leftward motion)              coherence            (rigthward motion)')
            ylabel('reaction time (ms)')
            title(sprintf('Reaction time data for session %d.', j), 'FontSize', 14);
            
        elseif task == 'd'
            xlabel('(CW rotation)                  disparity              (CCW rotation)')
            ylabel('reaction time ms')
            
            
        end %plot label
        
        
        title(sprintf('Reaction time data for session %d.', j), 'FontSize', 14);
        legend('correct choices', 'false choices')
        
%         figure 
%         hist(data(j).data(:,3));
%         title(sprintf('Reaction time histogram for session %d.', j), 'FontSize', 14);
        
    end %loop through sessions
    
    
else % combined data
    
    
    % transform rts as described above 
   rt_norm = log(data(:,3));
    
   sd = std(rt_norm);
   
   idx_rt_keep = rt_norm < (2.5 * sd);
   data_new = data(idx_rt_keep,:); 
   
   % select range over which we fit PMFs 
    if task == 'd'
        xrange = min(min(data_new(:,4)))-0.01:0.01:max(max(data_new(:,4)))+0.01; % range over which we evaluate fitted PMF
    elseif task == 'c'
        xrange = min(min(data_new(:,4)))-10:0.01:max(max(data_new(:,4)))+10;
        
    end
    
    % get data info we need for plot PMFs and median rts 
    [percent_right,n,confidence_interval,rt,coherences] = PMF_data(data_new);
   
    % code same as above 
    if Gauss
        
        
        
        % how to calculate initial guesses?
        x0 = [1,1]; % initial guesses for mean(x0(1)) and sd (x0(2))
        
        options = optimset(); % set options to default value - for fminsearch function
        
        [PMF_fit] = fminsearch(@cum_Gauss_loglike,x0,options,data_new);
        
        param.b0 = PMF_fit(1);
        param.b1 = PMF_fit(2);
        
        
        [p] = cum_Gauss_PMF(xrange,param,0);
        
        
        
    elseif logs
        
        
%         Y = [percent_right .* n ./100, n];
        
        Y = data_new(:,2);
        
%         X = [ones(length(coherences),1),coherences];
        
X = [ones(length(data_new(:,1)),1), data_new(:,4)];
        
        PMF_fit = glmfit(X,Y,'binomial','link','logit','constant','off');
        
        param.b0 = PMF_fit(1);
        param.b1 = PMF_fit(2);
        
        
        %
        [p] = logist_PMF(xrange,param,0);
        
        
        
        
    end %if Gauss combined data_new
    
    
    
    
    figure
    hold on
    plot(xrange,p.*100,'k-');
    errorbar(coherences,percent_right,confidence_interval,'ko')
    hold off
    ylim([0 100]);
    
    
    %plot label
    if task == 'c'
        xlabel('(leftward motion)              coherence            (rigthward motion)')
        ylabel('percentage rightward motion')
        title(sprintf('Data for combined sessions Mean = %.2f; SD = %.2f.', PMF_fit(1), PMF_fit(2)), 'FontSize', 14);
        
    elseif task == 'd'
        xlabel('(CW rotation)                  disparity              (CCW rotation)')
        ylabel('percentage CCW rotation')
        
    end %plot label
    
    %title
    if Gauss
        title(sprintf('Data for combined sessions Mean = %.2f; SD = %.2f.', PMF_fit(1), PMF_fit(2)), 'FontSize', 14);
    else
        title(sprintf('Data for combined sessions b0 = %.2f; b1 = %.2f.', PMF_fit(1), PMF_fit(2)), 'FontSize', 14);
    end %title
    
    
    
    %Plot Rts
    
    figure % correct rts
    errorbar(coherences,rt(:,1),rt(:,2),rt(:,5),'ko-')
    hold on
    errorbar(coherences,rt(:,3),rt(:,4),rt(:,6),'rd-')
    hold off
    %plot labels depending on rdk or cylinder task
    if task == 'c'
        xlabel('(leftward motion)              coherence            (rigthward motion)')
        ylabel('reaction time (ms)')
        
        
    elseif task == 'd'
        xlabel('(CW rotation)                  disparity              (CCW rotation)')
        ylabel('reaction time ms')
        
        
    end %plot label
    title('Reaction time data  combined session', 'FontSize', 14);
    legend('correct choices', 'false choices')
end %if data = struct




end %function






