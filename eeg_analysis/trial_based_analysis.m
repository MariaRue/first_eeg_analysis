


scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
     


    
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32];

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    clear data_append data 
for i = 1:6
cfg = []; 
cfg.dataset = fullfile(EEGdatadir,sprintf('fdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i)); 
cfg.trialdef.eventtype = 'trigger'; 
cfg.trialdef.eventvalue = [11, 30,40,50,130, 140, 150]; 
cfg.trialdef.prestim = 3;
cfg.trialdef.poststim = 10; 
cfg = ft_definetrial(cfg); 

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

coherence = [30,40,50];

figure

for bl = 1:4
    cfg = [];
  
    cfg.trials = data_append.trialinfo(:,2) == bl;
    
%     data_without_blinks.trialinfo(Ind(bl) : Ind(bl+1),2) = ones(length(Ind(bl) : Ind(bl+1)),1) * str2double(stim{1}.S.block_ID_cells{bl});
    data_select = ft_selectdata(cfg,data_without_blinks); 
for i = 1 : 3 % sort for coherences 
   
    idx_coh = data_select.trialinfo == coherence(i) | data_select.trialinfo == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_coherence{i} = ft_selectdata(cfg,data_select); 
    cfg = [];
average_ERP{i} = ft_timelockanalysis(cfg,data_coherence{i});

end

cfg = [];
cfg.channel = {'CPZ'};
cfg.baseline = [-1 -2];
cfg.baselinetype = 'absolute';
cfg.layout = 'quickcap64.mat';

% i think this is wrong and I plotted these differently in the actual analysis 
subplot(2,2,bl)

ft_singleplotER(cfg,average_ERP{1}, average_ERP{2}, average_ERP{3});



legend('30','40','50')

switch bl 
    
    case 1 
        condition = 'ITIS INTS'; 
    case 2 
        condition = 'ITIS INTL'; 
    case 3
        condition = 'ITIL INTS';
    case 4
        condition = 'ITIL INTL';
end 

title(sprintf('Subject: %d Condition: %s', subID, condition))
xlabel('secs, timelocked to trial start at 0')
% ft_singleplotER(cfg,average_ERP{3});
end
savefig(fullfile(EEGdir,'preprocessed_EEG_dat',sprintf('sub%03.0f_locked_to_tiral_start.fig',subID)))

end