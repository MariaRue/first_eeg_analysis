     % establish source density analysis and repeat O'connell style analyse to
% check whether we get the same results

scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
%EEGdir = '/Volumes/LaCie/data/';

addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer'));

addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28, 42];
%%
for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load data data_pre source_data
    
    if exist(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked'])) ~= 2
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_wo_blinks']));
    data = data_load.data_without_blinks;
    
    cfg = [];
    cfg.channel = {'all','-RM', '-VEOG', '-HEOG', '-LM'};
    [data_pre] = ft_preprocessing(cfg, data);
    % [easy_cap_labels] = change_electrode_labels(data.label);
    
    %data.label = easy_cap_labels;
    
    
    cfg = [];
    cfg.elec = data.elec;
    cfg.degree = 14;
    
    
    source_data =  ft_scalpcurrentdensity(cfg, data_pre);
    save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked']),'source_data');
    end
end
%% put all subjs into one dataframe

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked']));
    data{sj} = data_load.source_data;
    % change labels to easycap 
    [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    data{sj}.label = easy_cap_labels; 
     
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    data_avg =  compute_average_for_single_participant(data{sj},1, [-2 -1],0);
    coh_30{sj} = data_avg{1};
    coh_40{sj} = data_avg{2};
    coh_50{sj} = data_avg{3};
    
    data_avg_con =  compute_average_for_single_participant(data{sj},0, [-2 -1],0);
    cond_1{sj} = data_avg_con{1};
    cond_2{sj} = data_avg_con{2};
    cond_3{sj} = data_avg_con{3};
    cond_4{sj} = data_avg_con{4};
    
end

% lets calculate grand average for each coherence level 
cfg = []; 
coh_avg{1} = ft_timelockgrandaverage(cfg,coh_30{:});
coh_avg{2} = ft_timelockgrandaverage(cfg,coh_40{:});
coh_avg{3} = ft_timelockgrandaverage(cfg,coh_50{:});

% and for each condition 
cond_avg{1} = ft_timelockgrandaverage(cfg,cond_1{:});
cond_avg{2} = ft_timelockgrandaverage(cfg,cond_2{:});
cond_avg{3} = ft_timelockgrandaverage(cfg,cond_3{:});
cond_avg{4} = ft_timelockgrandaverage(cfg,cond_4{:});

 
%% average across subjects and conditions for each coherence level - timelocked to trial start 

lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);

 cl = cbrewer('seq','Blues',12);   
 cl =  cl([6 10 12],:);
minlim = -lim(2);
maxlim = lim(2);

cfg = [];
%cfg.channel = {'CPz'};
cfg.layout = 'easycapM1.mat';
cfg.ylim = [minlim maxlim];
% cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = cl;
cfg.linewidth = 3; 
%%%%%%%%%%%%


ft_singleplotER(cfg,coh_avg{:});


% YLIM - zero in middle!!! and make legend work!!!
% it *should have zero in the middle.
% and most importantly it should be the same in any related series of plots
% (e.g., all coherence levels)
legend({'30%', '40%', '50%'},'FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
set(gca,'FontSize',18)


%%  
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);

    
  
minlim = -lim(2);
maxlim = lim(2);
for i = 1:3
    
    start_time = -1;
    for t = 1:8
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [minlim maxlim];
cfg.layout = 'easycapM1.mat';
subplot(3,8,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,coh_avg{i}); colorbar
    end
end

subplot(3,8,1)
title('30% coherence timelocked to trial start at 0')
subplot(3,8,2)
title('0 1')
subplot(3,8,3)
title('1 2sec')
subplot(3,8,4)
title('2 3sec')
subplot(3,8,5)
title('3 4sec')
subplot(3,8,6)
title('4 5sec')
subplot(3,8,7)
title('5 6sec')
subplot(3,8,8)
title('6 7sec')

subplot(3,8,9)
title('40% coherence')

subplot(3,8,17)
title('50% coherence')
%%
for i = 1:24
    subplot(3,8,i)
    tidyfig;
    
end 
%% average across subjects and coherence levels for each condition - timelocked to trial start 





lim = quantile(cond_avg{1}.avg(:),[0.1 0.9]);

cl = cbrewer('div','RdBu', 12);  
cl = cl([4 1 9 12],:);
minlim = -10e-05;
maxlim = 10e-05;

cfg = [];
%cfg.channel = {'CPz'};
cfg.layout = 'easycapM1.mat';
cfg.ylim = [minlim maxlim];
% cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = cl;
cfg.linewidth = 3; 
%%%%%%%%%%%%







ft_singleplotER(cfg,cond_avg{:});

legend({'ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL'},'FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
set(gca,'FontSize',18)

%% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(cond_avg{1}.avg(:),[0.1 0.9]);

minlim = -lim(2);
maxlim = lim(2);
for i = 1:4
    start_time = -1;
    for t = 1:8
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [minlim maxlim];
cfg.layout = 'easycapM1.mat';
subplot(4,8,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,cond_avg{i}); colorbar
    end
end

subplot(4,8,1)
title('ITIS INTS condition timelocked to trial start at 0')
subplot(4,8,2)
title('0 1sec')
subplot(4,8,3)
title('1 2sec')
subplot(4,8,4)
title('2 3sec')
subplot(4,8,5)
title('3 4sec')
subplot(4,8,6)
title('4 5sec')
subplot(4,8,7)
title('5 6sec')

subplot(4,8,9)
title('ITIS INTL')

subplot(4,8,17)
title('ITIL INTS')

subplot(4,8,25)
title('ITIL INTL')

%% 
for i = 1:32
    subplot(4,8,i)
    tidyfig;
    
end

%% LRP 
% grand average between left and right motion topoplot for button press 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked']));
    data{sj} = data_load.source_data;
    
        [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    data{sj}.label = easy_cap_labels; 
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
% trial indicator

right_trials = data{sj}.trialinfo(:,1) == 30| data{sj}.trialinfo(:,1) == 40 | data{sj}.trialinfo(:,1) == 50;
left_trials = data{sj}.trialinfo(:,1) == 130| data{sj}.trialinfo(:,1) == 140 | data{sj}.trialinfo(:,1) == 150;

cfg = []; 
cfg.trials = right_trials;  
right_data = ft_selectdata(cfg,data{sj});
cfg = []; 
right_timelock{sj} = ft_timelockanalysis(cfg,right_data);


cfg = []; 
cfg.trials = left_trials; 
left_data = ft_selectdata(cfg,data{sj});
cfg = []; 
left_timelock{sj} = ft_timelockanalysis(cfg,left_data);
    



% calculate left - right grand average 

grand_average{sj} = left_timelock{sj}; 

grand_average{sj}.avg =left_timelock{sj}.avg - right_timelock{sj}.avg; 




 
 
 

end


cfg  = []; 
right_avg = ft_timelockgrandaverage(cfg,right_timelock{:}); 
left_avg = ft_timelockgrandaverage(cfg,left_timelock{:}); 
grand_average_avg = ft_timelockgrandaverage(cfg,grand_average{:}); 


cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [0 4]; 
cfg.zlim = [-6e-5 6e-5];
 ft_topoplotER(cfg,grand_average_avg); colorbar
 
 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
               right_lrp = ft_lateralizedpotential(cfg, left_avg,right_avg);
               
               
%% now calculate lrp for different coherence levels 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked']));
    data{sj} = data_load.source_data;
        [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    data{sj}.label = easy_cap_labels;     
    
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    
    data_avg_coh =  compute_average_for_single_participant(data{sj},1, [-2 -1],1);
    
    
   right_30{sj} = data_avg_coh{1,1};   
   right_40{sj} = data_avg_coh{2,1};  
   right_50{sj} = data_avg_coh{3,1};  
   
   left_30{sj} = data_avg_coh{1,2};   
   left_40{sj} = data_avg_coh{2,2};  
   left_50{sj} = data_avg_coh{3,2};  
 
 
end

cfg  = []; 
right_30_avg = ft_timelockgrandaverage(cfg,right_30{:}); 
right_40_avg = ft_timelockgrandaverage(cfg,right_40{:}); 
right_50_avg = ft_timelockgrandaverage(cfg,right_50{:}); 

cfg  = []; 
left_30_avg = ft_timelockgrandaverage(cfg,left_30{:}); 
left_40_avg = ft_timelockgrandaverage(cfg,left_40{:}); 
left_50_avg = ft_timelockgrandaverage(cfg,left_50{:}); 

grand_average_avg_30 = left_30_avg;
grand_average_avg_40 = left_40_avg;
grand_average_avg_50 = left_50_avg;

grand_average_avg_30.avg = left_30_avg.avg - right_30_avg.avg; 
grand_average_avg_40.avg = left_40_avg.avg - right_40_avg.avg; 
grand_average_avg_50.avg = left_50_avg.avg - right_50_avg.avg; 






cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [0 4]; 
cfg.zlim = [-10e-5 10e-5];

subplot(1,3,1)
ft_topoplotER(cfg,grand_average_avg_30); colorbar
title('30% coherence locked to trial start from trial start to +4sec')
set(gca,'FontSize',25)

subplot(1,3,2)
ft_topoplotER(cfg,grand_average_avg_40); colorbar
title('40%')
set(gca,'FontSize',25)

subplot(1,3,3)
ft_topoplotER(cfg,grand_average_avg_50); colorbar
title('50%')
set(gca,'FontSize',25)


cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
right_lrp_30 = ft_lateralizedpotential(cfg, left_30_avg,right_30_avg);
right_lrp_40 = ft_lateralizedpotential(cfg, left_40_avg,right_40_avg);
right_lrp_50 = ft_lateralizedpotential(cfg, left_50_avg,right_50_avg);
                

right_lrp_30.label = right_lrp_30.plotlabel;
right_lrp_40.label = right_lrp_40.plotlabel;
right_lrp_50.label = right_lrp_50.plotlabel;

 cl = cbrewer('seq','Blues',12);   
 cl =  cl([6 10 12],:);

figure
cfg = []; 
cfg.graphcolor = cl;
cfg.linewidth = 3;
ft_singleplotER(cfg,right_lrp_30, right_lrp_40,right_lrp_50);
legend({'30%', '40%', '50%'},'FontSize',25)

title('LRP timelocked to trial start for different coherences')
xlabel('time (secs), trial start at 0')
set(gca,'FontSize',25)




%%  repeat for the different conditions 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked']));
    data{sj} = data_load.source_data;
        [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    data{sj}.label = easy_cap_labels; 
    
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    
    data_avg_coh =  compute_average_for_single_participant(data{sj},0, [-2 -1],1);
    
    
   right_1{sj} = data_avg_coh{1,1};   
   right_2{sj} = data_avg_coh{2,1};  
   right_3{sj} = data_avg_coh{3,1};  
   right_4{sj} = data_avg_coh{4,1}; 
   
   left_1{sj} = data_avg_coh{1,2};   
   left_2{sj} = data_avg_coh{2,2};  
   left_3{sj} = data_avg_coh{3,2};  
   left_4{sj} = data_avg_coh{4,2}; 
 
end

cfg  = []; 
right_1_avg = ft_timelockgrandaverage(cfg,right_1{:}); 
right_2_avg = ft_timelockgrandaverage(cfg,right_2{:}); 
right_3_avg = ft_timelockgrandaverage(cfg,right_3{:}); 
right_4_avg = ft_timelockgrandaverage(cfg,right_4{:}); 

cfg  = []; 
left_1_avg = ft_timelockgrandaverage(cfg,left_1{:}); 
left_2_avg = ft_timelockgrandaverage(cfg,left_2{:}); 
left_3_avg = ft_timelockgrandaverage(cfg,left_3{:}); 
left_4_avg = ft_timelockgrandaverage(cfg,left_4{:});

grand_average_avg_1 = left_1_avg;
grand_average_avg_2 = left_2_avg;
grand_average_avg_3 = left_3_avg;
grand_average_avg_4 = left_4_avg;

grand_average_avg_1.avg = left_1_avg.avg - right_1_avg.avg; 
grand_average_avg_2.avg = left_2_avg.avg - right_2_avg.avg; 
grand_average_avg_3.avg = left_3_avg.avg - right_3_avg.avg; 
grand_average_avg_4.avg = left_4_avg.avg - right_4_avg.avg; 




cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [0 4]; 
cfg.zlim = [-10e-5 10e-5];

subplot(1,4,1)
ft_topoplotER(cfg,grand_average_avg_1); colorbar
title('ITIS INTS trial start to +4sec')
set(gca,'FontSize',25)

subplot(1,4,2)
ft_topoplotER(cfg,grand_average_avg_2); colorbar
title('ITIS INTL')
set(gca,'FontSize',25)

subplot(1,4,3)
ft_topoplotER(cfg,grand_average_avg_3); colorbar
title('ITIL INTS')
set(gca,'FontSize',25)

subplot(1,4,4)
ft_topoplotER(cfg,grand_average_avg_4); colorbar
title('ITIL INTL')
set(gca,'FontSize',25)

 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
right_lrp_1 = ft_lateralizedpotential(cfg, left_1_avg,right_1_avg);
right_lrp_2 = ft_lateralizedpotential(cfg, left_2_avg,right_2_avg);
right_lrp_3 = ft_lateralizedpotential(cfg, left_3_avg,right_3_avg);
right_lrp_4 = ft_lateralizedpotential(cfg, left_4_avg,right_4_avg);

right_lrp_1.label = right_lrp_1.plotlabel;
right_lrp_2.label = right_lrp_2.plotlabel;
right_lrp_3.label = right_lrp_3.plotlabel;
right_lrp_4.label = right_lrp_4.plotlabel;

cl = cbrewer('div','RdBu', 12);  
cl = cl([1 4 9 12],:);

figure;
cfg = []; 
cfg.graphcolor = cl;
cfg.linewidth = 3;
ft_singleplotER(cfg,right_lrp_1, right_lrp_2,right_lrp_3,right_lrp_4);
legend({'ITIS INTS%', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL'},'FontSize',25)

title('LRP timelocked to trial start for different conditions')
xlabel('time (secs), trial start at 0')
set(gca,'FontSize',25)


 



