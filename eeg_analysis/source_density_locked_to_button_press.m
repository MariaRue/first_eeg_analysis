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
    
    %if exist(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_trial_start_locked'])) ~= 2
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
    %end
end
%% 
clear data
for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']));
    data{sj} = data_load.source_data;
    
    
      
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    data_avg =  compute_average_for_single_participant(data{sj},1, [-6 -5],0);
    coh_30{sj} = data_avg{1};
    coh_40{sj} = data_avg{2};
    coh_50{sj} = data_avg{3};
    
    data_avg_con =  compute_average_for_single_participant(data{sj},0, [-6 -5],0);
    cond_1{sj} = data_avg_con{1};
    cond_2{sj} = data_avg_con{2};
    cond_3{sj} = data_avg_con{3};
    cond_4{sj} = data_avg_con{4};
    
end

% lets calculate grand average for each coherence level 
cfg = []; 
coh_avg_button_press{1} = ft_timelockgrandaverage(cfg,coh_30{:});
coh_avg_button_press{2} = ft_timelockgrandaverage(cfg,coh_40{:});
coh_avg_button_press{3} = ft_timelockgrandaverage(cfg,coh_50{:});

% and for each condition 
cond_avg_button_press{1} = ft_timelockgrandaverage(cfg,cond_1{:});
cond_avg_button_press{2} = ft_timelockgrandaverage(cfg,cond_2{:});
cond_avg_button_press{3} = ft_timelockgrandaverage(cfg,cond_3{:});
cond_avg_button_press{4} = ft_timelockgrandaverage(cfg,cond_4{:});

% change labels to easycap 
[easy_cap_labels] = change_electrode_labels(coh_avg_button_press{1}.label);


coh_avg_button_press{1}.label  = easy_cap_labels;
coh_avg_button_press{2}.label  = easy_cap_labels;
coh_avg_button_press{3}.label  = easy_cap_labels;
cond_avg_button_press{1}.label = easy_cap_labels;
cond_avg_button_press{2}.label = easy_cap_labels;
cond_avg_button_press{3}.label = easy_cap_labels;
cond_avg_button_press{4}.label = easy_cap_labels;
% 


%% average across subjects and conditions for each coherence level - timelocked to button press

lim = quantile(coh_avg_button_press{3}.avg(:),[0.1 0.9]);

 cl = cbrewer('seq','Blues',12);   
 cl =  cl([6 10 12],:);
% minlim = -lim(2);
% maxlim = lim(2);

maxlim = 30e-5;
minlim = -maxlim;


cfg = [];
%cfg.channel = {'CPz'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
% cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = cl;
cfg.linewidth = 3; 
%%%%%%%%%%%%


ft_singleplotER(cfg,coh_avg_button_press{:});


% YLIM - zero in middle!!! and make legend work!!!
% it *should have zero in the middle.
% and most importantly it should be the same in any related series of plots
% (e.g., all coherence levels)
legend({'30%', '40%', '50%'},'FontSize',25)
title('Averaged ERP across Subjects and conditions for different coherence levels CP1 CP2 CPZ P1 PZ P2','FontSize',25)
xlabel('time (s) - button press at 0','FontSize',14)
set(gca,'FontSize',25)

%%
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(coh_avg_button_press{3}.avg(:),[0.1 0.9]);

minlim = -lim(2);
maxlim = lim(2);
for i = 1:3
    start_time = -4;
    for t = 1:7
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-22e-5 22e-5];
cfg.layout = 'easycapM1.mat';
subplot(3,7,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,coh_avg_button_press{i}); colorbar
    end
end

subplot(3,7,1)
title('30% coherence timelocked to button press')
set(gca,'FontSize',25)



subplot(3,7,2)
title('-3 -2sec')
set(gca,'FontSize',25)
subplot(3,7,3)
title('-2 -1sec')
set(gca,'FontSize',25)
subplot(3,7,4)
title('-1 0sec')
set(gca,'FontSize',25)
subplot(3,7,5)
title('0 1sec')
set(gca,'FontSize',25)
subplot(3,7,6)
title('1 2sec')
set(gca,'FontSize',25)
subplot(3,7,7)
title('2 3sec')
set(gca,'FontSize',25)

subplot(3,7,8)
title('40% coherence')
set(gca,'FontSize',25)
subplot(3,7,15)
title('50% coherence')
set(gca,'FontSize',25)
%%
for i = 1:21
    subplot(3,7,i)
    tidyfig;
    
end 
%% %% average across subjects and coherence levels for each condition - timelocked to button press
lim = quantile(cond_avg_button_press{1}.avg(:),[0.1 0.9]);

cl = cbrewer('div','RdBu', 12);  
cl = cl([1 4 9 12],:);
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







ft_singleplotER(cfg,cond_avg_button_press{:});

legend({'ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL'},'FontSize',25)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
set(gca,'FontSize',25)
%% % now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(cond_avg_button_press{1}.avg(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);
for i = 1:4
    start_time = -4;
    for t = 1:7
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-22e-5 22e-5];
cfg.layout = 'easycapM1.mat';
subplot(4,7,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,cond_avg_button_press{i}); colorbar
    end
end

subplot(4,7,1)
title('ITIS INTS condition timelocked to button press')
set(gca,'FontSize',25)
subplot(4,7,2)
title('-3 -2sec')
set(gca,'FontSize',25)
subplot(4,7,3)
title('-2 -1sec')
set(gca,'FontSize',25)
subplot(4,7,4)
title('-1 0sec')
set(gca,'FontSize',25)
subplot(4,7,5)
title('0 1sec')
set(gca,'FontSize',25)
subplot(4,7,6)

title('1 2sec')
set(gca,'FontSize',25)
subplot(4,7,7)
title('2 3sec')
set(gca,'FontSize',25)
subplot(4,7,8)
title('ITIS INTL')
set(gca,'FontSize',25)
subplot(4,7,15)
title('ITIL INTS')
set(gca,'FontSize',25)
subplot(4,7,22)
title('ITIL INTL')
set(gca,'FontSize',25)
%%
for i = 1:28
    subplot(4,7,i)
    tidyfig;
    
end 
%% LRP 
% grand average between left and right motion topoplot for button press 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']));
    data{sj} = data_load.source_data;
    
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
right_avg_button_press = ft_timelockgrandaverage(cfg,right_timelock{:}); 
left_avg_button_press = ft_timelockgrandaverage(cfg,left_timelock{:}); 
grand_average_avg_button_press = ft_timelockgrandaverage(cfg,grand_average{:}); 

[easy_cap_labels] = change_electrode_labels(grand_average_avg_button_press.label);

grand_average_avg_button_press.label = easy_cap_labels; 

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [-2 2]; 
cfg.zlim = [-10e-5 10e-5];
 ft_topoplotER(cfg,grand_average_avg_button_press); colorbar
 
 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
               right_lrp = ft_lateralizedpotential(cfg, left_avg_button_press,right_avg_button_press);
 
 
 %% 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']));
    data{sj} = data_load.source_data;
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    
    data_avg_coh =  compute_average_for_single_participant(data{sj},1, [-6 -5],1);
    
    
   right_30{sj} = data_avg_coh{1,1};   
   right_40{sj} = data_avg_coh{2,1};  
   right_50{sj} = data_avg_coh{3,1};  
   
   left_30{sj} = data_avg_coh{1,2};   
   left_40{sj} = data_avg_coh{2,2};  
   left_50{sj} = data_avg_coh{3,2};  
 
 
end
%%
cfg  = []; 
right_30_avg_button_press = ft_timelockgrandaverage(cfg,right_30{:}); 
right_40_avg_button_press = ft_timelockgrandaverage(cfg,right_40{:}); 
right_50_avg_button_press = ft_timelockgrandaverage(cfg,right_50{:}); 

cfg  = []; 
left_30_avg_button_press = ft_timelockgrandaverage(cfg,left_30{:}); 
left_40_avg_button_press = ft_timelockgrandaverage(cfg,left_40{:}); 
left_50_avg_button_press = ft_timelockgrandaverage(cfg,left_50{:}); 

grand_average_avg_30_button_press = left_30_avg_button_press;
grand_average_avg_40_button_press = left_40_avg_button_press;
grand_average_avg_50_button_press = left_50_avg_button_press;

grand_average_avg_30_button_press.avg = left_30_avg_button_press.avg - right_30_avg_button_press.avg; 
grand_average_avg_40_button_press.avg = left_40_avg_button_press.avg - right_50_avg_button_press.avg; 
grand_average_avg_50_button_press.avg = left_50_avg_button_press.avg - right_40_avg_button_press.avg; 

[easy_cap_labels] = change_electrode_labels(grand_average_avg_30_button_press.label);

grand_average_avg_30_button_press.label = easy_cap_labels; 
grand_average_avg_40_button_press.label = easy_cap_labels; 
grand_average_avg_50_button_press.label = easy_cap_labels; 

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [-2 -2]; 
cfg.zlim = [-15e-5 15e-5];

subplot(1,3,1)
ft_topoplotER(cfg,grand_average_avg_30_button_press); colorbar
title('LRP button press -2 to +2 30% coherence')
set(gca,'FontSize',25)
subplot(1,3,2)
ft_topoplotER(cfg,grand_average_avg_40_button_press); colorbar
title('40%')
set(gca,'FontSize',25)
subplot(1,3,3)
ft_topoplotER(cfg,grand_average_avg_50_button_press); colorbar
title('50%')
set(gca,'FontSize',25)
 cfg = []; 
 cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
right_lrp_30_button_press = ft_lateralizedpotential(cfg, left_30_avg_button_press,right_30_avg_button_press);
right_lrp_40_button_press = ft_lateralizedpotential(cfg, left_40_avg_button_press,right_40_avg_button_press);
right_lrp_50_button_press = ft_lateralizedpotential(cfg, left_50_avg_button_press,right_50_avg_button_press);
                

right_lrp_30_button_press.label = right_lrp_30_button_press.plotlabel;
right_lrp_40_button_press.label = right_lrp_40_button_press.plotlabel;
right_lrp_50_button_press.label = right_lrp_50_button_press.plotlabel;

 cl = cbrewer('seq','Blues',12);   
 cl =  cl([6 10 12],:);
figure

cfg = []; 
cfg.graphcolor = cl;
cfg.linewidth = 3;
figure
ft_singleplotER(cfg,right_lrp_30_button_press, right_lrp_40_button_press,right_lrp_50_button_press);
legend({'30%', '40%', '50%'},'FontSize',25)

title('LRP timelocked to button press for different coherences')
xlabel('time (secs), button press at 0')
set(gca,'FontSize',25)


%%  repeat for the different conditions 

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load 
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_source_density_button_press_locked']));
    data{sj} = data_load.source_data;
    
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
%% 
cfg  = []; 
right_1_avg_button_press = ft_timelockgrandaverage(cfg,right_1{:}); 
right_2_avg_button_press = ft_timelockgrandaverage(cfg,right_2{:}); 
right_3_avg_button_press = ft_timelockgrandaverage(cfg,right_3{:}); 
right_4_avg_button_press = ft_timelockgrandaverage(cfg,right_4{:}); 

cfg  = []; 
left_1_avg_button_press = ft_timelockgrandaverage(cfg,left_1{:}); 
left_2_avg_button_press = ft_timelockgrandaverage(cfg,left_2{:}); 
left_3_avg_button_press = ft_timelockgrandaverage(cfg,left_3{:}); 
left_4_avg_button_press = ft_timelockgrandaverage(cfg,left_4{:});

grand_average_avg_1_button_press = left_1_avg_button_press;
grand_average_avg_2_button_press = left_2_avg_button_press;
grand_average_avg_3_button_press = left_3_avg_button_press;
grand_average_avg_4_button_press = left_4_avg_button_press;

grand_average_avg_1_button_press.avg = left_1_avg_button_press.avg - right_1_avg_button_press.avg; 
grand_average_avg_2_button_press.avg = left_2_avg_button_press.avg - right_2_avg_button_press.avg; 
grand_average_avg_3_button_press.avg = left_3_avg_button_press.avg - right_3_avg_button_press.avg; 
grand_average_avg_4_button_press.avg = left_4_avg_button_press.avg - right_4_avg_button_press.avg; 


[easy_cap_labels] = change_electrode_labels(grand_average_avg_4_button_press.label);
grand_average_avg_1_button_press.label = easy_cap_labels; 
grand_average_avg_2_button_press.label = easy_cap_labels; 
grand_average_avg_3_button_press.label = easy_cap_labels; 
grand_average_avg_4_button_press.label = easy_cap_labels; 


cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = [-2 2]; 
cfg.zlim = [-15e-5 15e-5];

subplot(1,4,1)
ft_topoplotER(cfg,grand_average_avg_1_button_press); colorbar
title('LRP button press -2 to +2 ITIS INTS')
set(gca,'FontSize',25)
subplot(1,4,2)
ft_topoplotER(cfg,grand_average_avg_2_button_press); colorbar
title('ITIS INTL')
subplot(1,4,3)
ft_topoplotER(cfg,grand_average_avg_3_button_press); colorbar
title('ITIL INTS')

subplot(1,4,4)
ft_topoplotER(cfg,grand_average_avg_4_button_press); colorbar
title('ITIL INTL')

cfg = []; 
cfg.channelcmb = {'C3' 'C4'
                   'C1' 'C2'
                   'CP3' 'CP4'
                   'CP1' 'CP2'}; 
               
               
right_lrp_1_button_press = ft_lateralizedpotential(cfg, left_1_avg_button_press,right_1_avg_button_press);
right_lrp_2_button_press = ft_lateralizedpotential(cfg, left_2_avg_button_press,right_2_avg_button_press);
right_lrp_3_button_press = ft_lateralizedpotential(cfg, left_3_avg_button_press,right_3_avg_button_press);
right_lrp_4_button_press = ft_lateralizedpotential(cfg, left_4_avg_button_press,right_4_avg_button_press);

right_lrp_1_button_press.label = right_lrp_1_button_press.plotlabel;
right_lrp_2_button_press.label = right_lrp_2_button_press.plotlabel;
right_lrp_3_button_press.label = right_lrp_3_button_press.plotlabel;
right_lrp_4_button_press.label = right_lrp_4_button_press.plotlabel;

cl = cbrewer('div','RdBu', 12);  
cl = cl([1 4 9 12],:);


cfg = []; 
cfg.graphcolor = cl;
cfg.linewidth = 3;

figure;
ft_singleplotER(cfg,right_lrp_1_button_press, right_lrp_2_button_press,right_lrp_3_button_press,right_lrp_4_button_press);
legend({'ITIS INTS%', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL'},'FontSize',25)

title('LRP timelocked to button press for different conditions')
xlabel('time (secs), button press at 0')
set(gca,'FontSize',25)




