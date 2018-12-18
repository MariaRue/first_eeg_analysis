% addpath('/Users/maria/Documents/data/data.continous_rdk/behavioural_pilot/sub002/');
% 
% Data = load ('sub002_sess010_behav.mat');

addpath('/Users/maria/Desktop/task_for_windows/behaviour'); 
Data = load ('sub000_sess067_behav.mat');

%% 
Stimulus = Data.S;
response = Data.respMat;


for i = 1:4
    idx_correct = response{i}(:,7) == 1;
    idx_incorrect = response{i}(:,7) == 0;
    idx_early = response{i}(:,7) == 2;
    idx_missed = response{i}(:,7) == 3;
    
    frame_correct = response{i}(idx_correct,6);
    frame_incorrect = response{i}(idx_incorrect,6);
    frame_early = response{i}(idx_early,6);
    frame_missed = response{i}(idx_missed,6);
    
    
    if Stimulus.block_ID_cells{i} == '1'
        
        t = 'ITI short, INTE short';
    elseif Stimulus.block_ID_cells{i} == '2'
        t = 'ITI short, INTE long';
        
    elseif Stimulus.block_ID_cells{i} == '3'
        
        t = 'ITI long, INTE short';
        
    elseif Stimulus.block_ID_cells{i} == '4'
        t = 'ITI long, INTE long';
    end
    
    subplot(4,1,i)
    plot(Stimulus.coherence_frame{i})
    hold on
    plot(Stimulus.mean_coherence{i})
    if i==9
        l(1) = plot(frame_correct,ones(numel(frame_correct),1),'g.')
        l(2) = plot(frame_incorrect,ones(numel(frame_incorrect),1),'rx')
        l(3) = plot(frame_missed,ones(numel(frame_missed),1),'kd')
        l(4) = plot(frame_early,ones(numel(frame_early),1),'bo')
        
        ylim([-1.5 1.5]);
        
        title(t)
        
        legend(l,{'correct','incorrect','missed','early'})
    else
        
        plot(frame_correct,ones(numel(frame_correct),1),'g.')
        plot(frame_incorrect,ones(numel(frame_incorrect),1),'rx')
        plot(frame_missed,ones(numel(frame_missed),1),'kd')
        plot(frame_early,ones(numel(frame_early),1),'bo')
        
        ylim([-1.5 1.5]);
        title(t)
    end
    hold off
    
    num_incoh_frames = sum(Stimulus.mean_coherence{i} == 0);

    ratio_early(i) = sum(idx_early)/num_incoh_frames .* 60 .* 60; %r resp per minute
    
    num_coh_frames = sum(Stimulus.mean_coherence_org{i} ~= 0);
    
    num_trials = max(max(Stimulus.blocks_shuffled{i}));
    
    ratio_missed(i) = (sum(idx_missed)/num_trials); % missed resp per minute
     
    idx_rts = ~isnan(response{i}(:,2));
   mean_rt(i) =  mean(response{i}(idx_rts,2));
end % loop through blocks
%% n

