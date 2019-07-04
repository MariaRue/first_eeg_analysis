


scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
     


    
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [28, 42];
%% get data timelocked to trial start 
for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    clear data_append data 
    
    if exist(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_wo_blinks.mat'])) ~= 2
for i = 1:6
cfg = []; 
cfg.dataset = fullfile(EEGdatadir,sprintf('fdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i)); 
cfg.trialdef.eventtype = 'trigger'; 
cfg.trialdef.eventvalue = [11, 30,40,50,130, 140, 150]; 
cfg.trialdef.prestim = 3;
cfg.trialdef.poststim = 10; 
cfg = ft_definetrial(cfg); 


  cfg.reref       = 'yes';
        cfg.channel     = 'all';
        cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
        cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
        cfg.detrend = 'yes';
        
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 5; % we played around with filtering to get rid of the
        % weird artefacts we found, usually also highpass filter with 0.1 and ord 3
        cfg.lpfiltord = 3;
data{i} = ft_preprocessing(cfg); 




 fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
        bhv{i} = load(fname_behav);
        fname_stim = fullfile(STdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
        stim{i} = load(fname_stim);
   Ind = find(data{i}.trialinfo(:,1) == 11);    
   
   data{i}.trialinfo(:,3) = ones(length(data{i}.trialinfo(:,1)),1) * i; 
        
 for bl = 1:4
    cfg = [];
    if bl <= 3
    cfg.trials = [Ind(bl) : Ind(bl+1)];
    data{i}.trialinfo(Ind(bl) : Ind(bl+1),2) = ones(length(Ind(bl) : Ind(bl+1)),1) * str2double(stim{i}.S.block_ID_cells{bl});
    else 
        cfg.trials = Ind(bl) : length(data{i}.trialinfo); 
        data{i}.trialinfo(Ind(bl) : length(data{i}.trialinfo),2) = ones(length(Ind(bl) : length(data{i}.trialinfo)),1) * str2double(stim{i}.S.block_ID_cells{bl});
    end
    
 end 
end
 cfg = []; 
data_append =  ft_appenddata(cfg,data{:});
   save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_append']),'data_append');
        
cd (scriptdir)

% remove eyeblinks 
data_without_blinks = data_append; 
for tr = 1:length(data_append.trial)
clear X
clear betas
clear predYblink 
clear Y 
 

X(:,1) =  data_append.trial{tr}(63,:); 
X(:,1) = X(:,1) - mean(X(:,1)); 
X(:,2) = ones(size(X(:,1)));

for i = 1:61

    Y(i,:) = data_append.trial{tr}(i,:); 
   
    betas(i,:) = glmfit(X,Y(i,:)','normal','constant','off'); 
    
end 

predYblink = betas(:,1)*X(:,1)'; 
%imagesc(Y - predYblink); 

data_without_blinks.trial{tr}(1:61,:) = Y - predYblink; 
end

 save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_wo_blinks']),'data_without_blinks');
    end
end 

%% put all subjs into one dataframe

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load
    
    
     data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trial_start_locked_wo_blinks']));
     data{sj} = data_load.data_without_blinks; 
%      coh_30_idx = data{sj}.trialinfo(:,1) == 30 | data{sj}.trialinfo(:,1) == 130;
%       coh_40_idx = data{sj}.trialinfo(:,1) == 40 | data{sj}.trialinfo(:,1) == 140;
%        coh_50_idx = data{sj}.trialinfo(:,1) == 50 | data{sj}.trialinfo(:,1) == 150;
%      
%        
%        cfg = []; 
%        cfg.trials = coh_30_idx; 
%        cfg.channel = 'CPZ'; 
%        cfg.avgoverchan = 'yes';
%        cfg.avgoverrpt = 'yes';
%      coh_30{sj} = ft_selectdata(cfg,data{sj}); 
%      
%             cfg = []; 
%        cfg.trials = coh_40_idx; 
%        cfg.channel = 'CPZ'; 
%      cfg.avgoverrpt = 'yes';
%      cfg.avgoverchan = 'yes';
%      coh_40{sj} = ft_selectdata(cfg,data{sj});  
%      
%             cfg = []; 
%        cfg.trials = coh_50_idx; 
%        cfg.channel = 'CPZ'; 
%      cfg.avgoverchan = 'yes';
%      cfg.avgoverrpt = 'yes';
%      coh_50{sj} = ft_selectdata(cfg,data{sj}); 
%      
%      
%         chan = 40; 
%    time = [2 3.5]; 
%   
%    timesl_coh = find( coh_30{sj}.time{1} >= time(1) &  coh_30{sj}.time{1} <= time(2));
%  
%     values_coh_30(sj)  = mean(coh_30{sj}.trial{1}(timesl_coh));
%   values_coh_40(sj)  = mean(coh_40{sj}.trial{1}(timesl_coh));
%     values_coh_50(sj)  = mean(coh_50{sj}.trial{1}(timesl_coh));
%      
     num_trials = length(data{sj}.trial);
     data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
     
    
      
end 

% 
% M1 = [values_coh_30',values_coh_40', values_coh_50'];
% figure; plot(M1','o-'); xlim([0.5 3.5])
% xticks([1 2 3])
% xticklabels({'30%','40%','50%'})
% legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
%         'subj7', 'subj8', 'subj9', 'subj10', 'subj11', 'sub12'}, 'location','EastOutside');
%  p_coh = anova1(M1);
 

cfg = []; 
data_all_subj = ft_appenddata(cfg,data{:});
[easy_cap_labels] = change_electrode_labels(data_all_subj.label);
data_all_subj.label = easy_cap_labels; 

%save(fullfile(EEGdir,'preprocessed_EEG_dat','all_subjs_trial_start'),'data_all_subj');
%% average across subjects and conditions for each coherence level - timelocked to trial start 
coherence = [30,40,50];

figure
for i = 1 : 3 % sort for coherences 
   
    idx_coh = data_all_subj.trialinfo(:,1) == coherence(i) | data_all_subj.trialinfo(:,1) == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_coherence{i} = ft_selectdata(cfg,data_all_subj); 
    
    cfg = []; 
    data_coherence{i} = ft_timelockanalysis(cfg, data_coherence{i});
    

    cfg = [];



average_ERP_timelock{i} = ft_timelockanalysis(cfg,data_coherence{i});

cfg = []; 
cfg.baseline = [-1 -2];
cfg.baselinetype = 'absolute';
average_ERP{i} = ft_timelockbaseline(cfg,average_ERP_timelock{i});
end



cfg = []; 
cfg.layout = 'easycapM1';

ft_singleplotER(cfg,average_ERP{1}, average_ERP{2}, average_ERP{3});
legend('30%', '40%', '50%','FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
%%  
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(average_ERP{1}.trial(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);
for i = 1:3
    start_time = -1;
    for t = 1:8
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-1.5 1.5];
cfg.layout = 'quickcap64.mat';
subplot(3,8,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,average_ERP{i}); colorbar
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

figure
for bl = 1 : 4 % sort for coherences 
   
    idx_coh = data_all_subj.trialinfo(:,2) == bl ;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_block{bl} = ft_selectdata(cfg,data_all_subj); 

    cfg = [];
cfg.channel = {'CPZ'};
cfg.baseline = [-1 -2];
cfg.baselinetype = 'absolute';

average_ERP{bl} = ft_timelockbaseline(cfg,data_block{bl});

end



cfg = []; 

ft_singleplotER(cfg,average_ERP{1}, average_ERP{2}, average_ERP{3}, average_ERP{4});

legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL','FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)

%%
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(average_ERP{1}.trial(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);
for i = 1:4
    start_time = -1;
    for t = 1:8
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-1.5 1.5];
cfg.layout = 'quickcap64.mat';
subplot(4,8,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,average_ERP{i}); colorbar
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

%% try plotting 