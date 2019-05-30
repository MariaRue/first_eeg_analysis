function rt_stair_distribution(data)
% empty vectors for combined data plot
rt_stair2_correct_all = []; % previous trial won right, current trial correct choice 
rt_stair3_correct_all = []; % previous trial won left, current trial correct choice 
rt_correct_all = [];
% rt_right_wrong_all = []; % previous trial won right, current trial incorrect choice
% rt_left_wrong_all = []; % previous trial won left, current trial incorrect choice 

if isstruct(data)
num_sessions = length(data); %number of sessions

for i = 1:num_sessions % loop through sessions
    mean_rt = [];
    se = [];
    
%     
    if length(data(i).data(:,1)) > 1000 
       data(i).data = data(i).data(1:1000,:);
     
    end 
% %  
    
    % remove outlier rts - transform to normal distrib with log, then take
    % 2.5sd as cut-off for outliers 
    
    rt_norm = log(data(i).data(:,3));
    
   sd = std(rt_norm);
   
   idx_rt_keep = rt_norm < (2.5 * sd);
%     
%     figure
%     hist(data(i).data(:,3));
%     title('untransformed')
%     
%     figure
%     hist(rt_norm);
%     title('transformed')

% figure 
% hist(data(i).data(idx_rt_keep,3));
% title('rt keep')


    
    data(i).data = data(i).data(idx_rt_keep,:);
    
     
    num_trials = length(data(i).data(:,1)); % number of trials in session
    
idx_stair2 = [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(2:num_trials,5) == 1)]; % last 2 trials won
idx_stair3 = [0; 0; (data(i).data(1:num_trials-2,5)== 1 & data(i).data(3:num_trials,5) == 1)];  % last 3 trials won 
idx_correct = [data(i).data(:,5)==1];
    
    %get rts for correct choices on current trial 
    rt_stair2_correct = data(i).data(logical(idx_stair2),3:4);
    rt_stair3_correct = data(i).data(logical(idx_stair3),3:4);
    rt_correct = data(i).data(logical(idx_correct),3:4);
 
   
    coherences = [];
    coherences = unique(data(i).data(:,4));
     
    
rt_stair2_correct_all = [rt_stair2_correct_all;rt_stair2_correct];
rt_stair3_correct_all = [rt_stair3_correct_all; rt_stair3_correct];
rt_correct_all = [rt_correct_all; rt_correct];


    
    %calculate mean rts and sds for incorrect and correct choices for
   % current trial]
   
   k = 1;
    for k = 1:length(coherences)
        
        mean_rt(k,1) = mean(rt_stair2_correct(rt_stair2_correct(:,2)==coherences(k),1)); %  mean rts for rightward correct choices
        mean_rt(k,2) = mean(rt_stair3_correct(rt_stair3_correct(:,2)==coherences(k),1)); %  mean rts for leftward correct choices
        mean_rt(k,3) = mean(rt_correct(rt_correct(:,2)==coherences(k),1)); %  mean rts for leftward correct choices
       
        
        se(k,1) = std(rt_stair2_correct(rt_stair2_correct(:,2)==coherences(k),1))/sqrt(length(rt_stair2_correct(rt_stair2_correct(:,2)==coherences(k),1)));
        se(k,2) = std(rt_stair3_correct(rt_stair3_correct(:,2)==coherences(k),1))/sqrt(length(rt_stair3_correct(rt_stair3_correct(:,2)==coherences(k),1))); %  ses for leftward correct choices
        se(k,3) = std(rt_correct(rt_correct(:,2)==coherences(k),1))/sqrt(length(rt_correct(rt_correct(:,2)==coherences(k),1))); %  ses for leftward correct choices
       
%         if k > 11
%             keyboard;
%         end
        
        
    end %loop through coherences
    
    
    mean_rt(isnan(mean_rt)) = 0;
    se(isnan(se)) = 0;
    
    figure
    hold on
    
    errorbar(coherences,mean_rt(:,1),se(:,1),'k')
    errorbar(coherences,mean_rt(:,2),se(:,2),'r')
    errorbar(coherences,mean_rt(:,3),se(:,3),'b')
    
    hold off
    
    legend('stair 2', 'stair 3','correct')
    title(sprintf('Mean rts with se for session = %.2f', i));
    xlabel('coherence')
    ylabel('rt')
    
    
  
    
    
end % loop through sessions



end % going through single sessions 

%plot combined data 

k = 1;
coherences = [];
coherences = unique(rt_stair2_correct_all(:,2));
    for k = 1:length(coherences)
        
        mean_rt_all(k,1) = mean(rt_stair2_correct_all(rt_stair2_correct_all(:,2)==coherences(k),1)); %  mean rts for rightward correct choices
        mean_rt_all(k,2) = mean(rt_stair3_correct_all(rt_stair3_correct_all(:,2)==coherences(k),1)); %  mean rts for leftward correct choices
        mean_rt_all(k,3) = mean(rt_correct_all(rt_correct_all(:,2)==coherences(k),1));
        
        se_all(k,1) = std(rt_stair2_correct_all(rt_stair2_correct_all(:,2)==coherences(k),1))/sqrt(length(rt_stair2_correct_all(rt_stair2_correct_all(:,2)==coherences(k),1)));
        se_all(k,2) = std(rt_stair3_correct_all(rt_stair3_correct_all(:,2)==coherences(k),1))/sqrt(length(rt_stair3_correct_all(rt_stair3_correct_all(:,2)==coherences(k),1))); %  ses for leftward correct choices
        se_all(k,3) = std(rt_correct_all(rt_correct_all(:,2)==coherences(k),1))/sqrt(length(rt_correct_all(rt_correct_all(:,2)==coherences(k),1))); %  ses for leftward correct choices
       
        
    end %loop through coherences 
    
    
     
        
    mean_rt(isnan(mean_rt)) = 0;
    se(isnan(se)) = 0;
    
    figure
    hold on
    
    errorbar(coherences,mean_rt_all(:,1),se_all(:,1),'k')
    errorbar(coherences,mean_rt_all(:,2),se_all(:,2),'r')
    errorbar(coherences,mean_rt_all(:,3),se_all(:,3),'b')
    
    hold off
    legend('stair2', 'stair3', 'correct')
    
    title('mean rts and se for combined data')
    xlabel('coherence')
    ylabel('rt')

    
    



end