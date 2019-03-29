% addpath('/Users/maria/Documents/data/data.continous_rdk/behavioural_pilot/sub002/');
% 
% Data = load ('sub002_sess010_behav.mat');

% BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub005/behaviour/'; 

BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/EEG/sub007/behaviour/'; 
Stimdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/EEG/sub007/stim/'; 
% BHVdatadir = '/Users/maria/Documents/data/data.continuous_rdk/data/training/sub010/behaviour/';
% Data = load ('sub004_sess001_behav.mat');


%% load in behavioural data
subID = 7; 
nSess = 6; 
for i = 1:nSess
   
    fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
    fname_sti = fullfile(Stimdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
    bhv{i} = load(fname_behav);
    stim{i} = load(fname_sti); 
end
%% 
% Stimulus = bhv{3}.S;
sess = 2; 
 response = bhv{sess}.respMat;



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
for i = 1:6
    
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
idx_3sec_INTEL_ITIS = poolRts_INTEL_ITIS <= 3.5; 
idx_3sec_INTEL_ITIL =  poolRts_INTEL_ITIL <= 3.5; 

poolRts_INTEL_ITIS = poolRts_INTEL_ITIS(idx_3sec_INTEL_ITIS);
 poolRts_INTEL_ITIL =  poolRts_INTEL_ITIL(idx_3sec_INTEL_ITIL);

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






