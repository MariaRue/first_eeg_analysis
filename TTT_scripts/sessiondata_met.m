function [session_data, cohlist] = sessiondata_met(trial_num,task)

% returns rt,choice, trial outcome, stimulus condition of session in matrix 

%Input:
%           trial_num = total number of trials in session - load
%                       sessdesc.mat file from a session - contains total trial number
%           task      = task id, 'd' for cylinder and 'c' for rdk Shadlen
%         

% This function is used in TTT_training_analysis_eye.m to extract the
% behavioural data from met structures and to save them in one big
% structure

%To run this function on its own one needs to navigate to
%the sessions trial folder eg.: /Users/dpag0617/Documents/data/data.nhp/
%TTT/behav.eye/M134.Turtle/20170518/M134.346.5.2AFC.rdk_motion.reactim/trials


%Output
% session_data = matrix with 5 columns

%Column     Information

%  1:       Trial condition - only relevant for tasks with a social or specific
%           reward condition - see TTT_training_analysis.m from Nela 
%           (left blank, was only inserted to make data matrix consistent 
%           across experiments (touch screen set up and eye tracking set up))

%  2:       choice -  0 = leftwards choice, 1 = rightwards choice  

%  3:       Reaction time (rt)

%  4:       Coherence/Disparity 

%  5:       choice - 0 = failed choice, 1 = correct choice 



% cohlist = list with disparity/coherence values 
  
% Maria Ruesseler, University of Oxford, May 2017


%load msig and parameter files if already in trial folder and analysing on
%single session 
for i = 1:trial_num
    
    istr = num2str(i); suff = ['_',istr,'.mat'];
    
    %parameter file
    td_in = load([istr,'/param',suff]); td(i) = td_in.td;
    
    %msig file
    msig(i)=load([istr,'/metsigs',suff]);
    
    %hit file
    hit(i) = load([istr,'/hitregion',suff]);
    
end % load msig and param file for a session


%counter for broken, ignored and aborted trials
broken_trial = 0;
%counter for matrix with trialinfo - coherence, rt, choice


trial = zeros(trial_num,5);
hitregion_vec = zeros(trial_num,1);
for n = 1:trial_num % looping through each trial of session and taking rt, choice and coherence, if trial has been completed
    
    % getting cargo that specifies whether trial has been completed or was
    % aborted, ignored or broken
    
    crg = msig(n).crg(msig(n).sig==3);
    
    %specifying whether trial has been completed or wheter it has not
    if crg >= 3
        broken_trial = broken_trial + 1;
        
        trial(n,:) = [NaN, NaN, NaN, NaN, NaN];
        hitregion_vec(n) = NaN;
    else
        % in case trial has been completed
        
        % get coherence/disparity 
        if task == 'c'
         trial(n,4) = td(n).stimlink(4).vpar.coherence * 100;
         
        elseif task == 'd'
        trial(n,4) = td(n).stimlink(2).vpar.disparity;
        
        
        else 
            fprintf('\n\Error: cannot determine task condition.n');
            
        end % get coherence/disparity 
        
        
        
        % get reaction time
       
        % possibly a different way to define the stop signal. Here I used
        % the time stamp at the time the programme identifies the eyes at
        % the target, but this incorporates the time the monkey actually
        % makes the saccade, which adds noise to the decision time. It
        % might be more accurate to go from this stop signal backwards and
        % time and look when the eye signal location changes for the first
        % time after presentation time of stimulus e.g. cylinder? 
            
            i = 7 == msig(n).sig  & (msig(n).crg == 4 | msig(n).crg == 5); % get trial stop signal when programme detects eyes at target location
            %l = 6 == msig(n).sig & msig(n).crg == 5; % get trial start signal when reaction time state starts 
            l = 6 == msig(n).sig & msig(n).crg == 4; % get trial start signal when stimulus appears 
            
      if sum(i) > 1 % Metbug = sometimes two mstop signals present in that case we always want the second stopsignal
          
          idx_t = find(i == 1);
          i(idx_t(1)) = 0;
          
            
      end
            %subtracting timestamps to get reaction time
            trial(n,3) = msig(n).tim(i) - msig(n).tim(l); %reaction time
  
        
        % get choice
        if msig(n).crg(i) == 4 % crg for failed choice
            trial(n,5) = 0;
            
        else
            trial(n,5) = 1; % crg for correct choice
            
        end % get choice 
        
        
        %correct target location 0 is on the left 1 on the right needed to calculate PMF data point for 0
        %coherence in which performance is shown as percentage
        %of rightward choices per coherence/disparity - based on Nelas
        %convention in training_analysis script
        
        % get the hitregion which indicates which target was correct 
        if task == 'c' % depending on task (cylinder or rdk) hitregion cell 4 or 5 have to be selected. hi is index
            hi = 4;
        elseif task == 'd'
            hi = 5;
        end
        
        
        % determining correct choice by looking at the location of the
        % hitregion from the MET software for each trial 
        if hit(n).hitregion{hi}(1) < 0
            hitregion_vec(n) = 0; % correct choice on the left 
        else
            hitregion_vec(n) = 1; % correct choice on the right 
        end % hitregion 
        
    end % if crg > 3
    
end % looping through each trial in a session



% figure out whether subjects choice was to the right or to the left 

% choices to the right 
rc = hitregion_vec == 1 & trial(:,5) == 1; % hitregion is to the right and monkeys choice was correct
lf = hitregion_vec == 0 & trial(:,5) == 0; % hitregion was left but monkeys choice was wrong so he made a choice to the right

%choices to the left
lc = hitregion_vec == 0 & trial(:,5)==1; % hitretion was on the left and monkey made a correct choice, i.e. to the left
rf = hitregion_vec == 1 & trial(:,5)==0; % hitregion was on the right and monkey made a wrong choice i.e. he chose the left target

% fill in choices in matrix column 
trial(rc|lf,2) = 1; % choices to the right 
trial(lc|rf,2) = 0; % choices to the left 

% remove NaN rows 
keep = isnan(trial)== 0;
session_data = trial(keep(:,1),:);

% get list of coherences 
tab_coherences = tabulate(trial(:,4));
cohlist = tab_coherences(:,1); % list of coherences, disparities that were present in session 



end % function



