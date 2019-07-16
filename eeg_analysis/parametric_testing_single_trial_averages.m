scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');




addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 53, 28, 42];

%% put all subjs into one dataframe

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_wo_blinks']));
    data{sj} = data_load.data_without_blinks;
    
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    for i = 1:6
        
        idx_coh  = data{sj}.trialinfo(:,3) == i;
        
        cfg = [];
        cfg.trials = idx_coh;
        data_sess{i} = ft_selectdata(cfg,data{sj});
        
        % divide each session for coherence
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,1) == 30 | data_sess{i}.trialinfo(:,1) == 130
        coh_30{i} = ft_selectdata(cfg, data_sess{i});
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        coh_30_time{i} = ft_timelockanalysis(cfg,coh_30{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        coh_30_time{i} = ft_timelockbaseline(cfg,coh_30_time{i});
        
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,1) == 40 | data_sess{i}.trialinfo(:,1) == 140
        coh_40{i} = ft_selectdata(cfg, data_sess{i});
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        coh_40_time{i} = ft_timelockanalysis(cfg,coh_40{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        coh_40_time{i} = ft_timelockbaseline(cfg,coh_40_time{i});
        
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,1) == 50 | data_sess{i}.trialinfo(:,1) == 150
        coh_50{i} = ft_selectdata(cfg, data_sess{i});
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        coh_50_time{i} = ft_timelockanalysis(cfg,coh_50{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        coh_50_time{i} = ft_timelockbaseline(cfg,coh_50_time{i});
        
        
        cfg = [];
        cfg.method = 'within';
        
        coh_30_avg_sj{sj} = ft_timelockgrandaverage(cfg, coh_30_time{:});
        coh_40_avg_sj{sj} = ft_timelockgrandaverage(cfg, coh_40_time{:});
        coh_50_avg_sj{sj} = ft_timelockgrandaverage(cfg, coh_50_time{:});
        
        %         cfg = [];
        %         cfg.baseline = [-2 -1];
        %         coh_30_avg_base_sj{sj} = ft_timelockbaseline(cfg,coh_30_avg_sj{sj});
        %         coh_40_avg_base_sj{sj} = ft_timelockbaseline(cfg,coh_40_avg_sj{sj});
        %         coh_50_avg_base_sj{sj} = ft_timelockbaseline(cfg,coh_50_avg_sj{sj});
        %         % divide each session for blocks
        
        
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,2) == 1;
        block_1{i} = ft_selectdata(cfg, data_sess{i});
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        block_1_time{i} = ft_timelockanalysis(cfg,block_1{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        block_1_time{i} = ft_timelockbaseline(cfg,block_1_time{i});
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,2) == 2;
        block_2{i} = ft_selectdata(cfg, data_sess{i});
        
        
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        block_2_time{i} = ft_timelockanalysis(cfg,block_2{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        block_2_time{i} = ft_timelockbaseline(cfg,block_2_time{i});
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,2) == 3;
        block_3{i} = ft_selectdata(cfg, data_sess{i});
        
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        block_3_time{i} = ft_timelockanalysis(cfg,block_3{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        block_3_time{i} = ft_timelockbaseline(cfg,block_3_time{i});
        
        cfg = [];
        cfg.trials = data_sess{i}.trialinfo(:,2) == 4;
        block_4{i} = ft_selectdata(cfg, data_sess{i});
        
        cfg = [];
        cfg.channel = 'all';
        cfg.baseline = [-1 -2];
        cfg.baselinetype = 'absolute';
        cfg.layout = 'quickcap64.mat';
        
        block_4_time{i} = ft_timelockanalysis(cfg,block_4{i});
        
        cfg = [];
        cfg.baseline = [-2 -1];
        block_4_time{i} = ft_timelockbaseline(cfg,block_4_time{i});
        
        
        cfg = [];
        cfg.method = 'within';
        
        block_1_avg_sj{sj} = ft_timelockgrandaverage(cfg, block_1_time{:});
        block_2_avg_sj{sj} = ft_timelockgrandaverage(cfg, block_2_time{:});
        block_3_avg_sj{sj} = ft_timelockgrandaverage(cfg, block_3_time{:});
        block_4_avg_sj{sj} = ft_timelockgrandaverage(cfg, block_4_time{:});
        %
        %
        %                         cfg = [];
        %         cfg.baseline = [-2 -1];
        %         block_1_avg_base_sj{sj} = ft_timelockbaseline(cfg,block_1_avg_sj{sj});
        %         block_2_avg_base_sj{sj} = ft_timelockbaseline(cfg,block_2_avg_sj{sj});
        %         block_3_avg_base_sj{sj} = ft_timelockbaseline(cfg,block_3_avg_sj{sj});
        %         block_4_avg_base_sj{sj} = ft_timelockbaseline(cfg,block_4_avg_sj{sj});
        %         figure (2)
        %         subplot(2,4,sj)
        %         cfg = [];
        %         cfg.channel = 'CPZ';
        %
        %         ft_singleplotER(cfg, block_1_avg_sj{sj},block_2_avg_sj{sj},block_3_avg_sj{sj},block_4_avg_sj{sj})
        %         if sj == 8
        %
        %             legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL','location','EastOutside')
        %         end
    end
    
    
    %
    figure
    cfg = [];
    cfg.channel = 'CPZ';
    ft_singleplotER(cfg,  coh_30_avg_sj{sj},  coh_40_avg_sj{sj},  coh_50_avg_sj{sj})
    
    if sj == 8
        
        legend('30', '40', '50', 'location','EastOutside')
    end
    
    
    
    
    %     figure
    %
    %
    %                 cfg = [];
    %                 cfg.channel = 'CPZ';
    %
    %                 ft_singleplotER(cfg, block_1_avg_sj{sj},block_2_avg_sj{sj},block_3_avg_sj{sj},block_4_avg_sj{sj})
    %                 if sj == 12
    %
    %                     legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL','location','EastOutside')
    %                 end
    %
    
end



%% % calculate and plot grandaverage across participants

% coherence
cfg = [];
cfg.parameter = 'avg';

coh_30_avg_all_sj = ft_timelockgrandaverage(cfg, coh_30_time{:});
coh_40_avg_all_sj = ft_timelockgrandaverage(cfg, coh_40_time{:});
coh_50_avg_all_sj = ft_timelockgrandaverage(cfg, coh_50_time{:});

%         cfg = [];
%         cfg.baseline = [-2 -1];
%         coh_30_avg_all_base_sj = ft_timelockbaseline(cfg,coh_30_avg_all_sj);
%         coh_40_avg_all_base_sj = ft_timelockbaseline(cfg,coh_40_avg_all_sj);
%         coh_50_avg_all_base_sj = ft_timelockbaseline(cfg,coh_50_avg_all_sj);
%

figure
cfg = [];
cfg.channel = 'CPZ';
ft_singleplotER(cfg,  coh_30_avg_all_sj,  coh_40_avg_all_sj,  coh_50_avg_all_sj)
legend('30%','40%', '50%','location','EastOutside')

%
    % condition
        cfg = [];
        cfg.parameter = 'avg';

         block_1_avg_all_sj = ft_timelockgrandaverage(cfg, block_1_time{:});
         block_2_avg_all_sj = ft_timelockgrandaverage(cfg, block_2_time{:});
         block_3_avg_all_sj = ft_timelockgrandaverage(cfg, block_3_time{:});
         block_4_avg_all_sj = ft_timelockgrandaverage(cfg, block_4_time{:});

%                  cfg = [];
%         cfg.baseline = [-2 -1];
%          block_1_avg_all_base_sj =  ft_timelockbaseline(cfg,block_1_avg_all_sj);
%          block_2_avg_all_base_sj =  ft_timelockbaseline(cfg,block_2_avg_all_sj);
%          block_3_avg_all_base_sj =  ft_timelockbaseline(cfg,block_3_avg_all_sj);
%          block_4_avg_all_base_sj =  ft_timelockbaseline(cfg,block_4_avg_all_sj);

                 figure
    cfg = [];
    cfg.channel = 'CPZ';
    ft_singleplotER(cfg,  block_1_avg_all_sj,  block_2_avg_all_sj,  block_3_avg_all_sj, block_4_avg_all_sj)
    
    legend('ITIS INTS','ITIS INTL', 'ITIL INTS', 'ITIL INTL','location','EastOutside')

%% plotting average values for time in which we expect difference in signal between conditions

chan = 40;
time = [2 3.5];

timesl_coh_30 = find( coh_30_avg_all_sj.time >= time(1) &  coh_30_avg_all_sj.time <= time(2));
timesl_coh_40 = find( coh_30_avg_all_sj.time >= time(1) &  coh_30_avg_all_sj.time <= time(2));
timesl_coh_50 = find( coh_30_avg_all_sj.time >= time(1) &  coh_30_avg_all_sj.time <= time(2));

% select the individual subject data from the time points and calculate the mean
for sj = 1:12
    values_coh_30(sj)  = mean(coh_30_avg_sj{sj}.avg(chan,timesl_coh_30));
    values_coh_40(sj)  = mean(coh_40_avg_sj{sj}.avg(chan,timesl_coh_40));
    values_coh_50(sj)  = mean(coh_50_avg_sj{sj}.avg(chan,timesl_coh_50));
end


M1 = [values_coh_30',values_coh_40', values_coh_50'];
figure; plot(M1','o-'); xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'30%','40%','50%'})
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
    'subj7', 'subj8', 'subj9', 'subj10', 'subj11', 'sub12'}, 'location','EastOutside');
p_coh = anova1(M1);

% condition

chan = 40;
time = [2 3.5];

timesl_cond_1 = find( block_1_avg_all_sj.time >= time(1) &  block_1_avg_all_sj.time <= time(2));
timesl_cond_2 = find( block_2_avg_all_sj.time >= time(1) &  block_2_avg_all_sj.time <= time(2));
timesl_cond_3 = find( block_3_avg_all_sj.time >= time(1) &  block_3_avg_all_sj.time <= time(2));
timesl_cond_4 = find( block_4_avg_all_sj.time >= time(1) &  block_4_avg_all_sj.time <= time(2));
% select the individual subject data from the time points and calculate the mean
for sj = 1:12
    values_cond_1(sj)  = mean(block_1_avg_sj{sj}.avg(chan,timesl_cond_1));
    values_cond_2(sj)  = mean(block_2_avg_sj{sj}.avg(chan,timesl_cond_2));
    values_cond_3(sj)  = mean(block_3_avg_sj{sj}.avg(chan,timesl_cond_3));
    values_cond_4(sj)  = mean(block_4_avg_sj{sj}.avg(chan,timesl_cond_4));
end


M2 = [values_cond_1',values_cond_2', values_cond_3',values_cond_4'];
figure; plot(M2','o-'); xlim([0.5 4.5])
xticks([1 2 3 4])
xticklabels({'ITIS INTS','ITIS INTL','ITIL INTS', 'ITIL INTL'})
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
    'subj7', 'subj8','subj9', 'subj10', 'subj11','sub12'}, 'location','EastOutside');

Y = [M2(:,1); M2(:,2); M2(:,3); M2(:,4)];
% Itis(1) and ITIL(2) factor

g1(1:24) = 1;

g1(25:48) = 2;


% INTL (4) INTS (3)
g2(1:12) = 3;
g2(13:24) = 4;
g2(25:36) = 3;
g2(37:48) = 4;

p_cond = anovan(Y,{g1, g2});
%%


