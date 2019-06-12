
scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');

addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32];
%%
for sj = 1:length(subj_list)
    clear data_append data data_without_blinks data
    
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    
    for i = 1:6
        cfg = [];
        cfg.dataset = fullfile(EEGdatadir,sprintf('fdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = [11,201,205];
        cfg.trialdef.prestim = 7;
        cfg.trialdef.poststim = 5;
        
        
        
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
        % test whether num of button presses is the same for eeg and behav
        behav_trig = [];
        
        for bl = 1:4
            behav_trig = [behav_trig; bhv{i}.B.trigger_vals{bl}];
        end
        if sum(behav_trig == 201 | behav_trig == 205) ~= (length(data{i}.trial)-sum(data{i}.trialinfo(:,1)==11))
            
            keyboard;
            
        else
            % at some point account for wrong trials (but they barely happend and
            % are ignored here (button press in wrong motion direction)
            for bl = 1:4
                coherence = [];
                % find frame of button press
                frames_behav =  find(bhv{i}.B.trigger_vals{bl} == 201 | bhv{i}.B.trigger_vals{bl} == 205);
                
                
                for f = 1:length(frames_behav)
                    
                    % find frame in mean coherence
                    mean_coherence = stim{i}.S.mean_coherence_org{bl}(frames_behav(f));
                    t = 1;
                    while mean_coherence == 0
                        
                        mean_coherence = stim{i}.S.mean_coherence_org{bl}(frames_behav(f)-t);
                        t = t+1;
                    end
                    
                    if mean_coherence > 0
                        mean_coherence = (mean_coherence * 100) + 100;
                        
                    else
                        mean_coherence = mean_coherence * -100;
                    end
                    coherence = [coherence; mean_coherence];
                    
                    
                    
                end
                
                if bl < 4
                    data{i}.trialinfo(Ind(bl)+1 : Ind(bl+1)-1,3) = coherence;
                else
                    data{i}.trialinfo(Ind(bl)+1 : end,3) = coherence;
                    
                end
                
            end
            
            
        end
        
        
        
    end
    
    cfg = [];
    data_append =  ft_appenddata(cfg,data{:});
    save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_button_press_locked_append']),'data_append');
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
    
    save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_button_press_locked_wo_blinks']),'data_without_blinks');
    
end

%% put all subjs into one dataframe

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load
    
    
    data_load = load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_button_press_locked_wo_blinks']));
    data{sj} = data_load.data_without_blinks;
    num_trials = length(data{sj}.trial);
    data{sj}.trialinfo(: ,4) = ones(num_trials,1) * subID;
    
    
    
end

cfg = [];
data_all_subj = ft_appenddata(cfg,data{:});
save(fullfile(EEGdir,'preprocessed_EEG_dat','all_subjs_button_press'),'data_all_subj');

%% average across subjects and conditions for each coherence level - timelocked to trial start
coherence = [30,40,50];

figure
for i = 1 : 3 % sort for coherences
    
    idx_coh = data_all_subj.trialinfo(:,3) == coherence(i) | data_all_subj.trialinfo(:,3) == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh;
    data_coherence{i} = ft_selectdata(cfg,data_all_subj);
    
    cfg = [];
    cfg.channel = {'CPZ'};
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
    cfg.layout = 'quickcap64.mat';
    average_ERP{i} = ft_timelockbaseline(cfg,data_coherence{i});
    
end

ft_singleplotER(cfg,average_ERP{1}, average_ERP{2}, average_ERP{3});
legend('30%', '40%', '50%','FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
%%
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(average_ERP{1}.trial(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);
for i = 1:3
    start_time = -4;
    for t = 1:7
cfg = [];
cfg.xlim = [start_time start_time + 1];
start_time = start_time + 1; 
cfg.zlim = [-5 4];
cfg.layout = 'quickcap64.mat';
subplot(3,7,sb_idx)
sb_idx = sb_idx + 1; 
 ft_topoplotER(cfg,average_ERP{i}); colorbar
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
%% %% average across subjects and coherence levels for each condition - timelocked to trial start

figure
for bl = 1 : 4 % sort for coherences
    
    idx_coh = data_all_subj.trialinfo(:,2) == bl ;
    
    cfg = [];
    cfg.trials = idx_coh;
    data_block{bl} = ft_selectdata(cfg,data_all_subj);
    
    cfg = [];
    cfg.channel = {'CPZ'};
    cfg.baseline = [-6 -5];
    cfg.baselinetype = 'absolute';
    cfg.layout = 'quickcap64.mat';
    average_ERP{bl} = ft_timelockbaseline(cfg,data_block{bl});
    
end





ft_singleplotER(cfg,average_ERP{1}, average_ERP{2}, average_ERP{3}, average_ERP{4});

legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL','FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
%%
% now plot topoplot 
figure; 
sb_idx = 1; 

lim = quantile(average_ERP{1}.trial(:),[0.1 0.9]);

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
 ft_topoplotER(cfg,average_ERP{i}); colorbar
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