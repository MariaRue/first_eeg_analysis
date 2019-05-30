function rt_distributions(data)
% plot rts with regard to won t-1 left or right 

% empty vectors for combined data plot
rt_right_correct_all = []; % previous trial won right, current trial correct choice 
rt_left_correct_all = []; % previous trial won left, current trial correct choice 
rt_right_wrong_all = []; % previous trial won right, current trial incorrect choice
rt_left_wrong_all = []; % previous trial won left, current trial incorrect choice 

if isstruct(data)
num_sessions = length(data); %number of sessions

for i = 1:num_sessions % loop through sessions
    mean_rt = [];
    se = [];
    
    
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
    
    winright =  [0; (data(i).data(1:num_trials-1,5)== 1 & data(i).data(1:num_trials-1,2) == 1)]; % previous trial won right
    winleft = [0; (data(i).data(1:end-1,5)== 1 & data(i).data(1:end-1,2) == 0)]; % previous trial won left%   
    
    %correct choice on current trial
    idx_right_correct = winright == 1 & data(i).data(:,5) == 1;
    idx_left_correct = winleft == 1 & data(i).data(:,5) == 1;
    
    %incorrect choice on current trial 
    idx_right_wrong = winright == 1 & data(i).data(:,5) == 0;
    idx_left_wrong = winleft == 1 & data(i).data(:,5) == 0; 
    
    %get rts for correct choices on current trial 
    rt_right_correct = data(i).data(idx_right_correct,3:4);
    rt_left_correct = data(i).data(idx_left_correct,3:4);
 
    %get rts for incorrect choices on current trial 
    rt_right_wrong = data(i).data(idx_right_wrong,3:4);
    rt_left_wrong = data(i).data(idx_left_wrong,3:4);
    
    
    coherences = unique(data(i).data(:,4));
    
    
     % get rts for correct choices on current trial 
    rt_right_correct = rt_right_correct(rt_right_correct(:,1) <= 1,:);
    rt_left_correct = rt_left_correct(rt_left_correct(:,1) <= 1,:);
 
    % get rts for incorrect choices on current trial 
    rt_right_wrong = rt_right_wrong(rt_right_wrong(:,1) <= 1,:);
    rt_left_wrong = rt_left_wrong(rt_left_wrong(:,1) <= 1,:);
    
    
 
rt_right_correct_all = [rt_right_correct_all;rt_right_correct];
rt_left_correct_all = [rt_left_correct_all; rt_left_correct]; 

rt_right_wrong_all = [rt_right_wrong_all; rt_right_wrong]; 
rt_left_wrong_all = [rt_left_wrong_all; rt_left_wrong];  
    
    
    %calculate mean rts and sds for incorrect and correct choices for
   % current trial
    for k = 1:length(coherences)
        
        mean_rt(k,1) = mean(rt_right_correct(rt_right_correct(:,2)==coherences(k),1)).* 1000; %  mean rts for rightward correct choices
        mean_rt(k,2) = mean(rt_left_correct(rt_left_correct(:,2)==coherences(k),1)).* 1000; %  mean rts for leftward correct choices
        
        mean_rt(k,3) = mean(rt_right_wrong(rt_right_wrong(:,2)==coherences(k),1)).* 1000; % same as above but for failed choices
        mean_rt(k,4) = mean(rt_left_wrong(rt_left_wrong(:,2)==coherences(k),1)).* 1000;
        
        
      

        se(k,1) = std(rt_right_correct(rt_right_correct(:,2)==coherences(k),1))/sqrt(length(rt_right_correct(rt_right_correct(:,2)==coherences(k),1)));
        se(k,2) = std(rt_left_correct(rt_left_correct(:,2)==coherences(k),1))/sqrt(length(rt_left_correct(rt_left_correct(:,2)==coherences(k),1))); %  ses for leftward correct choices
        
        se(k,3) = std(rt_right_wrong(rt_right_wrong(:,2)==coherences(k),1))/sqrt(length(rt_right_wrong(rt_right_wrong(:,2)==coherences(k),1))); % same as above but for failed choices
        se(k,4) = std(rt_left_wrong(rt_left_wrong(:,2)==coherences(k),1))/sqrt(length(rt_left_wrong(rt_left_wrong(:,2)==coherences(k),1)));
% %         
%         if k > 11
%             keyboard;
%         end
        
        
    end %loop through coherences
    
    
    mean_rt(isnan(mean_rt)) = 0;
    se(isnan(se)) = 0;
    
%     figure
%     hold on
%     
%     errorbar(coherences,mean_rt(:,1),se(:,1),'k')
%     errorbar(coherences,mean_rt(:,2),se(:,2),'r')
%     errorbar(coherences,mean_rt(:,3),se(:,3),'b')
%     errorbar(coherences,mean_rt(:,4),se(:,4),'g')
%     
%     hold off
%     
%     legend('right correct', 'left correct', 'right wrong', 'left wrong' )
%     title(sprintf('Mean rts with se for session = %.2f', i));
%     xlabel('coherence')
%     ylabel('rt')
    
    
  
    
    
end % loop through sessions



end % going through single sessions 

%plot combined data 


coherences = unique(rt_right_correct_all(:,2));
    for k = 1:length(coherences)
        
        mean_rt_all(k,1) = median(rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1)); %  mean rts for rightward correct choices
        mean_rt_all(k,2) = median(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1)); %  mean rts for leftward correct choices
        
        mean_rt_all(k,3) = median(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1)); % same as above but for failed choices
        mean_rt_all(k,4) = median(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1));
% 
%         se_all(k,1) = prctile(rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1),25);
%         se_all(k,2) = prctile(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1),25); %  ses for leftward correct choices
%         
%         se_all(k,3) = prctile(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1),25); % same as above but for failed choices
%         se_all(k,4) = prctile(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1),25);
%          
        se_all(k,1) = mean_rt_all(k,1) - quantile([rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1)],0.25);
        se_all(k,3) = mean_rt_all(k,2) - quantile(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1),0.25); %  ses for leftward correct choices
%         
        se_all(k,5) = mean_rt_all(k,3)- quantile(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1),0.25); % same as above but for failed choices
        se_all(k,7) = mean_rt_all(k,4) - quantile(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1),0.25);
        
        se_all(k,2) = quantile(rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1),0.75) - mean_rt_all(k,1);
        se_all(k,4) = quantile(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1),0.75) - mean_rt_all(k,2); %  ses for leftward correct choices
%         
        se_all(k,6) = quantile(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1),0.75)-mean_rt_all(k,3); % same as above but for failed choices
        se_all(k,8) = quantile(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1),0.75)-mean_rt_all(k,4);

        
%         se_all(k,1) = std(rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1))/sqrt(length(rt_right_correct_all(rt_right_correct_all(:,2)==coherences(k),1)));
%         se_all(k,2) = std(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1))/sqrt(length(rt_left_correct_all(rt_left_correct_all(:,2)==coherences(k),1))); %  ses for leftward correct choices
%         
%         se_all(k,3) = std(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1))/sqrt(length(rt_right_wrong_all(rt_right_wrong_all(:,2)==coherences(k),1))); % same as above but for failed choices
%         se_all(k,4) = std(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1))/sqrt(length(rt_left_wrong_all(rt_left_wrong_all(:,2)==coherences(k),1)));
%         
%         
        
    end %loop through coherences 
    
    
     
        
    mean_rt(isnan(mean_rt)) = 0;
    se(isnan(se)) = 0;
    
    figure
    hold on
    
    errorbar(coherences,mean_rt_all(:,1),se_all(:,1),se_all(:,2),'k')
    errorbar(coherences,mean_rt_all(:,2),se_all(:,3),se_all(:,4),'r')
    hold off 
     legend('t-1 won right, t won', 't-1 won left, t won') 
    
     
     
     figure
    hold on 
    errorbar(coherences,mean_rt_all(:,3),se_all(:,5),se_all(:,6),'b')
    errorbar(coherences,mean_rt_all(:,4),se_all(:,7),se_all(:,8),'g')

%     plot(coherences,mean_rt_all(:,1),'k')
%    plot(coherences,mean_rt_all(:,2),'r')
%     plot(coherences,mean_rt_all(:,3),'b')
%     plot(coherences,mean_rt_all(:,4),'g')
%     
    hold off
   legend( 't-1 won right, t lost', 't-1 won left, t lost' )
    
    title('mean rts and se for combined data')
    xlabel('coherence')
    ylabel('rt')

    
    



end