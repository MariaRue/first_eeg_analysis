% establish source density analysis and repeat O'connell style analyse to
% check whether we get the same results

scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
%EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
EEGdir = '/Volumes/LaCie/data/';



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
  
%     coh_30_idx = data{sj}.trialinfo(:,1) == 30 | data{sj}.trialinfo(:,1) == 130;
%     coh_40_idx = data{sj}.trialinfo(:,1) == 40 | data{sj}.trialinfo(:,1) == 140;
%     coh_50_idx = data{sj}.trialinfo(:,1) == 50 | data{sj}.trialinfo(:,1) == 150;
%     
%     
%     cfg = [];
%     cfg.trials = coh_30_idx;
%     cfg.channel = 'CPZ';
%     cfg.avgoverchan = 'yes';
%     cfg.avgoverrpt = 'yes';
%     coh_30{sj} = ft_selectdata(cfg,data{sj});
%     
%     cfg = [];
%     cfg.trials = coh_40_idx;
%     cfg.channel = 'CPZ';
%     cfg.avgoverrpt = 'yes';
%     cfg.avgoverchan = 'yes';
%     coh_40{sj} = ft_selectdata(cfg,data{sj});
%     
%     cfg = [];
%     cfg.trials = coh_50_idx;
%     cfg.channel = 'CPZ';
%     cfg.avgoverchan = 'yes';
%     cfg.avgoverrpt = 'yes';
%     coh_50{sj} = ft_selectdata(cfg,data{sj});
%     
%     
%     chan = 40;
%     time = [2 3.5];
%     
%     timesl_coh = find( coh_30{sj}.time{1} >= time(1) &  coh_30{sj}.time{1} <= time(2));
%     
%     values_coh_30(sj)  = mean(coh_30{sj}.trial{1}(timesl_coh));
%     values_coh_40(sj)  = mean(coh_40{sj}.trial{1}(timesl_coh));
%     values_coh_50(sj)  = mean(coh_50{sj}.trial{1}(timesl_coh));
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    data_avg =  compute_average_for_single_participant(data{sj},1, [-2 -1]);
    coh_30{sj} = data_avg{1};
    coh_40{sj} = data_avg{2};
    coh_50{sj} = data_avg{3};
    
    data_avg_con =  compute_average_for_single_participant(data{sj},0, [-2 -1]);
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

% change labels to easycap 
[easy_cap_labels] = change_electrode_labels(coh_avg{1}.label);


coh_avg{1}.label  = easy_cap_labels;
coh_avg{2}.label  = easy_cap_labels;
coh_avg{3}.label  = easy_cap_labels;
cond_avg{1}.label = easy_cap_labels;
cond_avg{2}.label = easy_cap_labels;
cond_avg{3}.label = easy_cap_labels;
cond_avg{4}.label = easy_cap_labels;
% 
% M1 = [values_coh_30',values_coh_40', values_coh_50'];
% figure; plot(M1','o-'); xlim([0.5 3.5])
% xticks([1 2 3])
% xticklabels({'30%','40%','50%'})
% legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
%     'subj7', 'subj8', 'subj9', 'subj10', 'subj11', 'sub12'}, 'location','EastOutside');
% p_coh = anova1(M1);


% cfg = [];

% save(fullfile(EEGdir,'preprocessed_EEG_dat','all_subjs_trial_start'),'data_all_subj'); 
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
title('30% coherence timelocked to button press')
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




%%%%%%%%%%%%
% TOM SAYS %
%%%%%%%%%%%%
% when singleplotting, try to choose sensible colours that allow us to
% infer the 2 x 2 design (e.g., light blue, dark blue, light red, dark red)

% try something like this
% % get colours from colourbrewer
% cl = cbrewer('div','RdBu', 12);
% % select a dark red, a light red, a light blue, and a dark blue
% cl = cl([1 4 9 12],:);
% 
% % you might need to change the order of the rows of 'cl' around to it
% plots how you want... eg by doing: cl = cl([1 2 4 3], :);

% then when you plot use the following arguments
% cfg.graphcolor = cl; % calls nice colours
% cfg.linewidth = 3; % makes lines wider so they are easier to see from far
% away
% you can also use set(gca,'FontSize', 18) or some larger number, so the
% axis labels are larger

lim = quantile(cond_avg{1}.avg(:),[0.1 0.9]);

cl = cbrewer('div','RdBu', 12);  
% cl = cl([1 4 9 12],:);
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

lim = quantile(block_average_ERP{1}.avg(:),[0.1 0.9]);

minlim = lim(1);
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
 ft_topoplotER(cfg,block_average_ERP{i}); colorbar
    end
end

subplot(4,8,1)
title('ITIS INTS condition timelocked to button press')
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


% trial indicator 
right_trials = data_all_subj.trialinfo(:,1) == 30| data_all_subj.trialinfo(:,1) == 40 | data_all_subj.trialinfo(:,1) == 50;
left_trials =  data_all_subj.trialinfo(:,1) == 130| data_all_subj.trialinfo(:,1) == 140 | data_all_subj.trialinfo(:,1) == 150;

cfg = []; 
cfg.trials = right_trials;  
right_data = ft_selectdata(cfg,data_all_subj);
cfg = []; 
right_timelock = ft_timelockanalysis(cfg,right_data);

    
    cfg = [];
 
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
right_baseline = ft_timelockbaseline(cfg, right_timelock);

cfg = []; 
cfg.trials = left_trials; 
left_data = ft_selectdata(cfg,data_all_subj);
cfg = []; 
left_timelock = ft_timelockanalysis(cfg,left_data);
    
    cfg = [];
 
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
left_baseline = ft_timelockbaseline(cfg, left_timelock);


% calculate left - right grand average 

grand_average = left_baseline; 

grand_average.avg =left_baseline.avg - right_baseline.avg; 


cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [-1 1]; 

 ft_topoplotER(cfg,grand_average); colorbar
 
 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
right_lrp = ft_lateralizedpotential(cfg, left_baseline,right_baseline);
 
 
 


