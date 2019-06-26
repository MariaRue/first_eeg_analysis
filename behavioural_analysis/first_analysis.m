% addpath('/Users/maria/Documents/data/data.continous_rdk/behavioural_pilot/sub002/');
% 
% Data = load ('sub002_sess010_behav.mat');

% BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub005/behaviour/'; 
session = 'train';
switch session
    
    case 'EEG'
        BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/EEG/sub019/behaviour/';
        Stimdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/EEG/sub019/stim/';
        
    case 'train'
        BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/training/sub020/behaviour/';
        Stimdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/training/sub020/stim/';
end

% BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/training/sub010/behaviour/';
% Data = load ('sub004_sess001_behav.mat');
% BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/behaviour/sub000'
% Stimdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/stim/sub000'
%% load in behavioural data
 subID = 20; 
 nSess = 2;

sess = [10:11];
for i = 1:nSess
   s = sess(i);
    fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,s));
    fname_sti = fullfile(Stimdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,s));
    bhv{i} = load(fname_behav);
   stim{i} = load(fname_sti); 
end
%% 
% Stimulus = bhv{3}.S;
sess = 2; 
 response = bhv{sess}.respMat;


figure
for i = 1:4
    idx_correct = response{i}(:,7) == 1;
    idx_incorrect = response{i}(:,7) == 0;
    idx_early = response{i}(:,7) == 2;
    idx_missed = response{i}(:,7) == 3;
    

    
    frame_correct = response{i}(idx_correct,6);
    frame_incorrect = response{i}(idx_incorrect,6);
    frame_early = response{i}(idx_early,6);
    frame_missed = response{i}(idx_missed,6);
    
    
    if stim{sess}.S.block_ID_cells{i} == '1'
        
        t = 'ITI short, INTE short';
    elseif stim{sess}.S.block_ID_cells{i} == '2'
        t = 'ITI short, INTE long';
        
    elseif stim{sess}.S.block_ID_cells{i} == '3'
        
        t = 'ITI long, INTE short';
        
    elseif stim{sess}.S.block_ID_cells{i} == '4'
        t = 'ITI long, INTE long';
    end
     
    subplot(4,1,i)
    plot(bhv{sess}.B.coherence_frame{i})
    hold on
    plot(bhv{sess}.B.mean_coherence{i})
  
    
    if i==9
        l(1) = plot(frame_correct,ones(numel(frame_correct),1),'g.')
        l(2) = plot(frame_incorrect,ones(numel(frame_incorrect),1),'rx')
        l(3) = plot(frame_missed,ones(numel(frame_missed),1),'kd')
        l(4) = plot(frame_early,ones(numel(frame_early),1),'bo')
        
        ylim([-1.5 1.5]);
        
        title(t)
        
        legend(l,{'correct','incorrect','missed','early'})
    else
        
        plot(frame_correct,ones(numel(frame_correct),1),'g.','MarkerSize',15)
        plot(frame_incorrect,ones(numel(frame_incorrect),1),'rx')
        plot(frame_missed,ones(numel(frame_missed),1),'kd')
        plot(frame_early,ones(numel(frame_early),1),'bo')
        
        try
        if sum(idx_incorrect) == 0
             legend({'coherence', 'mean coherence', 'correct','missed','early'})
        else 
         legend({'coherence', 'mean coherence', 'correct','incorrect','missed','early'})
        end 
        catch
            %executed if error
        end
        
        ylim([-1.5 1.5]);
        title(t)
    end
    hold off
    
    num_incoh_frames = sum(bhv{sess}.B.mean_coherence{i} == 0);

    ratio_early(i) = sum(idx_early)/num_incoh_frames .* 60 .* 60; %r resp per minute
    
    num_coh_frames = sum(stim{sess}.S.mean_coherence_org{i} ~= 0);
    
    num_trials = max(max(stim{sess}.S.blocks_shuffled{i}));
    
    ratio_missed(i) = (sum(idx_missed)/num_trials); % missed resp per minute
     
    idx_rts = ~isnan(response{i}(:,2));
   mean_rt(i) =  mean(response{i}(idx_rts,2));
end % loop through blocks
%% look at RTs 

poolRts_INTEL_ITIS = []; 
poolRts_INTEL_ITIL = []; 
poolRts_INTES_ITIL = []; 
poolRts_INTES_ITIS = []; 
for i = 1:2
    
    clear response 
    response = bhv{i}.respMat; 
for b  = 1:4

    blockID = str2double(stim{i}.S.block_ID_cells{b}); 
switch blockID
    
    case 1 
        poolRts_INTES_ITIS = [poolRts_INTES_ITIS; response{b}(response{b}(:,7)==1,2)];
        
    case 2 
        poolRts_INTEL_ITIS =  [poolRts_INTEL_ITIS; response{b}(response{b}(:,7)==1,2)]; 
        
    case 3 
        poolRts_INTES_ITIL = [poolRts_INTES_ITIL; response{b}(response{b}(:,7)==1,2)]; 
    case 4 
       poolRts_INTEL_ITIL = [poolRts_INTEL_ITIL; response{b}(response{b}(:,7)==1,2)]; 
end 

end 
end


% only using trials with rths below 3.5 seconds for INTEL - to compare with
% INTES condition and see whether people adopt behaviour 

% idx_3sec_INTEL_ITIS = poolRts_INTEL_ITIS <= 3.5; 
% idx_3sec_INTEL_ITIL =  poolRts_INTEL_ITIL <= 3.5; 
% 
% poolRts_INTEL_ITIS = poolRts_INTEL_ITIS(idx_3sec_INTEL_ITIS);
%  poolRts_INTEL_ITIL =  poolRts_INTEL_ITIL(idx_3sec_INTEL_ITIL);

figure
subplot(1,2,1)
hold on
histogram(poolRts_INTEL_ITIL,10)
histogram(poolRts_INTES_ITIL,10)
hold off 
legend('INTEL ITIL','INTES ITIL')
title(sprintf('sub %d', subID))
xlabel('Rts (sec)')

subplot(1,2,2)
hold on
histogram(poolRts_INTEL_ITIS,10)
histogram(poolRts_INTES_ITIS,10)
hold off 
legend('INTEL ITIS','INTES ITIS')

%% look at Rts but now extend it to button presses made after the 3.5 cutoff in the short trial conditions 

% to get these rts loop through each trial that has a missed response and
% find the first button press after the missed response (make sure that it
% is an early button response) and calculate the rt time 

% another important note about this - should I disregard Rts below 0.5 secs
% and do the same for those trials as well? Should I check whether the
% button press response was correct in case of the early button press ones?
% 
poolRts_INTEL_ITIS = []; 
poolRts_INTEL_ITIL = []; 
poolRts_INTES_ITIL = []; 
poolRts_INTES_ITIS = []; 
figure;
for i = 1: 6
       
    clear response 
    response = bhv{i}.respMat; 
for bl  = 1:4
RT = [];
  
    idx_correct = response{bl}(:,7) == 1;
    idx_incorrect = response{bl}(:,7) == 0;
    idx_early = response{bl}(:,7) == 2;
    idx_missed = response{bl}(:,7) == 3;

    
    frame_correct = response{bl}(idx_correct,6);
    frame_incorrect = response{bl}(idx_incorrect,6);
    frame_early = response{bl}(idx_early,6);
    frame_missed = response{bl}(idx_missed,6);
    
    
    % find indices of trial starts 
    idx = bhv{1}.B.mean_coherence{bl}(2:end) ~= 0 & bhv{1}.B.mean_coherence{bl}(1:end-1) == 0; 
    
    idx_sot = find(idx); 
       blockID = str2double(stim{i}.S.block_ID_cells{bl}); 
    % loop through trials and identify the ones that were missed 
    if any(idx_missed)
    tr_id = []; 
    count_missed = 1; 
    for tr = 1:length(idx_sot)
        
     f_ID = frame_missed(count_missed); 
     
     diff_trs = diff([idx_sot(tr),f_ID]); 
     

   
   if blockID == 2 || blockID == 4
       tr_diff = [550 499];
       
   else
       tr_diff = [350 299];
       
   end 
     
     if  diff_trs > tr_diff(1)
     
        tr_id(tr) = 0; 
        
     elseif diff_trs < tr_diff(2)
         tr_id(tr) = 0;
        
     elseif diff_trs < tr_diff(1)
        tr_id(tr) = 1; 
        
        if length(frame_missed) > count_missed
        count_missed = count_missed + 1; 
        end
        
        
     % else 
     
        
        
    end 
    
    end
    

    
    % loop through trials that have been missed and find next early button
    % press 
    
    idx_missed_trial_id = find(tr_id); 
    count = 1; 
    for tr = 1:length(idx_missed_trial_id) 
        
     
     trial =   idx_sot(idx_missed_trial_id(tr)); 
     
     
    early_trials = find(frame_early > trial); 
    
    
    if any(early_trials)
    % now find early press that is closest to trial start (which is first
    % entry in early trials 
    
   % find the difference in frames 
 if idx_sot(idx_missed_trial_id(tr)) < idx_sot(end)
     
   if frame_early(early_trials(1)) < idx_sot(idx_missed_trial_id(tr)+1) 
       
   diff_frames = frame_early(early_trials(1))-trial; 
   RT(count) = diff_frames/100; 
   
   end 
   
 elseif  idx_sot(idx_missed_trial_id(tr)) == idx_sot(end)
     
   diff_frames = frame_early(early_trials(1))-trial; 
   RT(count) = diff_frames/100; 
   
     
   
 end 
 
  count = count + 1;
    end 
        
        
    end 
    
   RT(RT == 0) = [];
    end
    
    
   
switch blockID
    
    case 1 
        poolRts_INTES_ITIS = [poolRts_INTES_ITIS; response{bl}(response{bl}(:,7)==1,2);RT'];
        
    case 2 
        poolRts_INTEL_ITIS =  [poolRts_INTEL_ITIS; response{bl}(response{bl}(:,7)==1,2);RT']; 
        
    case 3 
        poolRts_INTES_ITIL = [poolRts_INTES_ITIL; response{bl}(response{bl}(:,7)==1,2);RT']; 
    case 4 
       poolRts_INTEL_ITIL = [poolRts_INTEL_ITIL; response{bl}(response{bl}(:,7)==1,2);RT']; 
end 

end 
end



subplot(1,2,1)
hold on
histogram(poolRts_INTEL_ITIL, 20, 'BinWidth',0.2,'BinLimits',[0 10])
histogram(poolRts_INTES_ITIL, 20, 'BinWidth',0.2,'BinLimits',[0 10])
hold off 
legend('INTEL ITIL','INTES ITIL')
title(sprintf('sub %d', subID))
xlabel('Rts (sec)')

subplot(1,2,2)
hold on
histogram(poolRts_INTEL_ITIS,20, 'BinWidth',0.2,'BinLimits',[0 10])
histogram(poolRts_INTES_ITIS,20, 'BinWidth',0.2, 'BinLimits',[0 10])
hold off 
legend('INTEL ITIS','INTES ITIS')
%%
    RT = [];
    idx_correct = response{bl}(:,7) == 1;
    idx_incorrect = response{bl}(:,7) == 0;
    idx_early = response{bl}(:,7) == 2;
    idx_missed = response{bl}(:,7) == 3;

    
    frame_correct = response{bl}(idx_correct,6);
    frame_incorrect = response{bl}(idx_incorrect,6);
    frame_early = response{bl}(idx_early,6);
    frame_missed = response{bl}(idx_missed,6);
    
    
    % find indices of trial starts 
    idx = bhv{1}.B.mean_coherence{bl}(2:end) ~= 0 & bhv{1}.B.mean_coherence{bl}(1:end-1) == 0; 
    
    idx_sot = find(idx); 
    
    % loop through trials and identify the ones that were missed 
    
    tr_id = []; 
    count_missed = 1; 
    for tr = 1:length(idx_sot)
        
     f_ID = frame_missed(count_missed); 
     
     diff_trs = diff([idx_sot(tr),f_ID]); 
     
  
     if  diff_trs > 350
     
        tr_id(tr) = 0; 
        
     elseif diff_trs < 299 
         tr_id(tr) = 0;
        
     elseif diff_trs < 350 
        tr_id(tr) = 1; 
        
        if length(frame_missed) > count_missed
        count_missed = count_missed + 1; 
        end
     else 
     
        
        
    end 
    
    end
    
    % loop through trials that have been missed and find next early button
    % press 
    
    idx_missed_trial_id = find(tr_id); 
    count = 1; 
    for tr = 1:length(idx_missed_trial_id) 
        
     
     trial =   idx_sot(idx_missed_trial_id(tr)); 
     
     
    early_trials = find(frame_early > trial); 
    
    
    if any(early_trials)
    % now find early press that is closest to trial start (which is first
    % entry in early trials 
    
   % find the difference in frames 
 if idx_sot(idx_missed_trial_id(tr)) < idx_sot(end)
     
   if frame_early(early_trials(1)) < idx_sot(idx_missed_trial_id(tr)+1) 
       
   diff_frames = frame_early(early_trials(1))-trial; 
   RT(count) = diff_frames/100; 
   
   end 
   
 elseif  idx_sot(idx_missed_trial_id(tr)) == idx_sot(end)
     
   diff_frames = frame_early(early_trials(1))-trial; 
   RT(count) = diff_frames/100; 
   
     
   
 end 
 
  count = count + 1;
    end 
        
        
    end 
    
   RT(RT == 0) = [];
   RT(RT > 8) = [];
   RT = [];
    
    blockID = str2double(stim{i}.S.block_ID_cells{bl}); 
switch blockID
    
    case 1 
        poolRts_INTES_ITIS = [poolRts_INTES_ITIS; response{bl}(response{bl}(:,7)==1,2);RT'];
        
    case 2 
        poolRts_INTEL_ITIS =  [poolRts_INTEL_ITIS; response{bl}(response{bl}(:,7)==1,2);RT']; 
        
    case 3 
        poolRts_INTES_ITIL = [poolRts_INTES_ITIL; response{bl}(response{bl}(:,7)==1,2);RT']; 
    case 4 
       poolRts_INTEL_ITIL = [poolRts_INTEL_ITIL; response{bl}(response{bl}(:,7)==1,2);RT']; 
end 

    
    

  
    
figure

subplot(1,2,1)
hold on
histogram(poolRts_INTEL_ITIL,'BinWidth',0.2)
histogram(poolRts_INTES_ITIL,'BinWidth',0.2)
hold off 
legend('INTEL ITIL','INTES ITIL')
title(sprintf('sub %d', subID))
xlabel('Rts (sec)')

subplot(1,2,2)
hold on
histogram(poolRts_INTEL_ITIS,20)
histogram(poolRts_INTES_ITIS,20)
hold off 
legend('INTEL ITIS','INTES ITIS')

%% 



