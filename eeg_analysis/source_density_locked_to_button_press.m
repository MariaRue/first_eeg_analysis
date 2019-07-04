% establish source density analysis and repeat O'connell style analyse to
% check whether we get the same results

scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
%EEGdir = '/Volumes/LaCie/data/';



addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28, 42];
%% source density analysis locked to button press 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load data data_pre source_data
    
    if exist(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked'])) ~= 2
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_button_press_locked_wo_blinks']));
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
    save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']),'source_data');
    end
end
%% put all subjs into one dataframe (button press)
clear data
for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']));
    data{sj} = data_load.source_data;
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    
    
end

cfg = [];
data_all_subj_button = ft_appenddata(cfg,data{:});
[easy_cap_labels] = change_electrode_labels(data_all_subj_button.label);
data_all_subj_button.label = easy_cap_labels; 
%save(fullfile(EEGdir,'preprocessed_EEG_dat','all_subjs_button_press'),'data_all_subj');

%% average across subjects and conditions for each coherence level - timelocked to button press
coherence = [30,40,50]; 

figure
for i = 1 : 3 % sort for coherences
    
    idx_coh = data_all_subj_button.trialinfo(:,3) == coherence(i) | data_all_subj_button.trialinfo(:,3) == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh;
    data_coherence{i} = ft_selectdata(cfg,data_all_subj_button);
    
    cfg = [];
%     cfg.channel = {'CPZ'};
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
    cfg.layout = 'quickcap64.mat';
    button_average_ERP_time{i} = ft_timelockanalysis(cfg,data_coherence{i});
    
    cfg = []; 
    
    button_average_ERP{i} = ft_timelockbaseline(cfg,button_average_ERP_time{i});
end

cfg = []; 
cfg.legend = {'30%', '40%', '50%'};
 cfg.layout = 'easycapM1.mat';
ft_singleplotER(cfg,button_average_ERP{1}, button_average_ERP{2}, button_average_ERP{3});
legend('30%', '40%', '50%')
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
%%
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(button_average_ERP_time{i}.avg(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);
for i = 1:3
    start_time = -4;
    for t = 1:7
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-5 4];
cfg.layout = 'easycapM1.mat';
subplot(3,7,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,button_average_ERP{i}); colorbar
    end
end

subplot(3,7,1)
title('30% coherence timelocked to button press')
subplot(3,7,2)
title('-3 -2sec')
subplot(3,7,3)
title('-2 -1sec')
subplot(3,7,4)
title('-1 0sec')
subplot(3,7,5)
title('0 1sec')
subplot(3,7,6)
title('1 2sec')
subplot(3,7,7)
title('2 3sec')

subplot(3,7,8)
title('40% coherence')

subplot(3,7,15)
title('50% coherence')
%%
for i = 1:21
    subplot(3,7,i)
    tidyfig;
    
end 
%% %% average across subjects and coherence levels for each condition - timelocked to button press
figure
clear average_ERP_time
clear average_ERP
for bl = 1 : 4 % sort for coherences
    
    idx_coh = data_all_subj_button.trialinfo(:,2) == bl ;
    
    cfg = [];
    cfg.trials = idx_coh;
    data_block{bl} = ft_selectdata(cfg,data_all_subj_button);
    
    cfg = [];
 
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
    cfg.layout = 'quickcap64.mat';
    block_button_average_ERP_time{bl} = ft_timelockanalysis(cfg,data_block{bl});
    block_button_average_ERP{bl} = ft_timelockbaseline(cfg,block_button_average_ERP_time{bl});
end



cfg = []; 
 cfg.layout = 'easycapM1.mat';
ft_singleplotER(cfg,block_button_average_ERP{1}, block_button_average_ERP{2}, block_button_average_ERP{3}, block_button_average_ERP{4});

legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL','FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
%% % now plot topoplot 
figure; 
sb_idx = 1; 

% lim = quantile(average_ERP_time{1}.avg(:),[0.1 0.9]);
lim = [-1.2516    0.6446];
minlim = lim(1);
maxlim = lim(2);
for i = 1:4
    start_time = -4;
    for t = 1:7
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-5 4];
cfg.layout = 'quickcap64.mat';
subplot(4,7,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,block_button_average_ERP{i}); colorbar
    end
end

subplot(4,7,1)
title('ITIS INTS condition timelocked to button press')
subplot(4,7,2)
title('-3 -2sec')
subplot(4,7,3)
title('-2 -1sec')
subplot(4,7,4)
title('-1 0sec')
subplot(4,7,5)
title('0 1sec')
subplot(4,7,6)
title('1 2sec')
subplot(4,7,7)
title('2 3sec')

subplot(4,7,8)
title('ITIS INTL')

subplot(4,7,15)
title('ITIL INTS')

subplot(4,7,22)
title('ITIL INTL')
%%
for i = 1:28
    subplot(4,7,i)
    tidyfig;
    
end 
%% LRP 
% grand average between left and right motion topoplot for button press 


% trial indicator 
right_trials = data_all_subj_button.trialinfo(:,3) == 30| data_all_subj_button.trialinfo(:,3) == 40 | data_all_subj_button.trialinfo(:,3) == 50;
left_trials =  data_all_subj_button.trialinfo(:,3) == 130| data_all_subj_button.trialinfo(:,3) == 140 | data_all_subj_button.trialinfo(:,3) == 150;

cfg = []; 
cfg.trials = right_trials;  
right_data = ft_selectdata(cfg,data_all_subj_button);
cfg = []; 
button_right_timelock = ft_timelockanalysis(cfg,right_data);

    
    cfg = [];
 
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
button_right_baseline = ft_timelockbaseline(cfg, button_right_timelock);

cfg = []; 
cfg.trials = left_trials; 
left_data = ft_selectdata(cfg,data_all_subj_button);
cfg = []; 
button_left_timelock = ft_timelockanalysis(cfg,left_data);
    
    cfg = [];
 
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
button_left_baseline = ft_timelockbaseline(cfg, button_left_timelock);


% calculate left - right grand average 

button_grand_average = button_left_baseline; 

button_grand_average.avg = button_left_baseline.avg - button_right_baseline.avg; 


cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [-1 1]; 

 ft_topoplotER(cfg,button_grand_average); colorbar
 
 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
 button_lrp = ft_lateralizedpotential(cfg, button_left_baseline, button_right_baseline);
 
 
 


