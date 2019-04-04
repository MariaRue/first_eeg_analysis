%% This script is a newer version of how prelim the preprocessing of the EEG data could work (13Dec 2018) 

% data will be read in as a continuous version, then eyeblinks will be
% removed by regressing the EOG channel out of the raw data. After we try
% to identify the electrical artifacts found in the EEG data seem to
% reflect big voltage jumps 

%% start fieldtrip
addpath('/Users/maria/Documents/Matlab/fieldtrip/')
ft_defaults;
%% 
addpath('/Users/maria/Documents/eeglab14_1_2b/')
eeglab

%% path to data
addpath('/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub003/EEG');
eegdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub003/EEG';

%%
 load('/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub003/behaviour/sub003_sess004_behav.mat')

%% load data JUST TO FIND OUT WHERE THE EPOCH BEGINS AND ENDS (in terms of samples)

cfg                        = [];
cfg.dataset                = fullfile(eegdatadir,'sub003_sess002_fil001.set');
cfg.trialdef.eventtype     = 'trigger';

cfg.trialdef.eventvalue    = [24 25 26 34 35 36 210]; % your event values
cfg.trialdef.prestim       = 1;  % before stimulation (sec), only use positive num

cfg.trialdef.poststim      = 2; % after stimulation (sec) , only use positive num

cfg                        = ft_definetrial(cfg);

% we don't need the data, just the trial info
% data   = ft_preprocessing(cfg);

% keep the trial info for later
cfg_short_trials = cfg;

%% cut out the entire epoch as one trial


end_of_block = find(cfg_short_trials.trl(:,4)==210);

for bl = 1 : length(end_of_block)
cfg = [];
cfg.dataset                 = fullfile(eegdatadir,'sub003_sess002_fil001.set');
% cfg.trialdef.eventtype      = 'trigger';
% first and last samples


if bl == 1

cfg.trl                     = [cfg_short_trials.trl(1) cfg_short_trials.trl(end_of_block(bl),2) 0];
else 
 
    cfg.trl                     = [cfg_short_trials.trl(end_of_block(bl-1)+1,1) cfg_short_trials.trl(end_of_block(bl),2) 0];
end 
% 
% filter 
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;

data_one_big_epoch{bl}          = ft_preprocessing(cfg);
end 
%% downsample 
% cfg = []; 
% cfg.resamplefs = 150;
% data_down = ft_resampledata(cfg,data_one_big_epoch);

data_down = data_one_big_epoch{1}; 


%% 

cfg = []; 
cfg.viewmode = 'vertical'; 
cfg.channel = 'EEG'; 

cfg = ft_databrowser(cfg, data_down); 
%% regress eyeblinks out 

clear X
clear betas
clear predYblink 
clear Y 
clear data_without_blinks 

X(:,1) =  data_down.trial{1}(63,:); 
X(:,1) = X(:,1) - mean(X(:,1)); 
X(:,2) = ones(size(X(:,1)));

for i = 1:61

    Y(i,:) = data_down.trial{1}(i,:); 
   
    betas(i,:) = glmfit(X,Y(i,:)','normal','constant','off'); 
    
end 

predYblink = betas(:,1)*X(:,1)'; 
imagesc(Y - predYblink); 

data_without_blinks = Y - predYblink; 

figure; 
subplot(2,1,1) 
imagesc(Y);
title('continuous data before eyeblink removal') 
ylabel('EEG channell')
xlabel('time') 
subplot(2,1,2)
imagesc(Y-predYblink);
title('continuous data after eyeblink removal')

figure; 
c = 1; 
for i = [2, 40]; 
subplot(4,1,c) 
plot(data_down.trial{1}(i,:)) 

ylabel('voltage') 
xlabel('time') 


c = c+1; 
subplot(4,1,c) 
plot(data_down.time{1},data_without_blinks(i,:)) 

c = c+1; 

end
%% remove noisy channels 

data_eye = data_down; 
data_eye.trial{1}(1:61,:) = data_without_blinks(); 


   cfg.keepchannel = 'yes'; 
   cfg.channel = 'EEG';

[data_clean] =  ft_rejectvisual(cfg,data_eye);


%%  remove EEG artifacts 
clear z_data; 
% z-score data 
for i = 1:length(data_clean.trial{1}(:,1)) 
   
    z_data(i,:) = (data_clean.trial{1}(i,:) - mean(data_clean.trial{1}(i,:)))/std(data_clean.trial{1}(i,:)); 
    
    
    
end 

figure 
subplot(2,1,1) 
plot(z_data(1,:)) 

subplot(2,1,2) 
plot(z_data(35,:)) 


%% replace artifacts with NaNs 
clear smoothed_abs_zdata
clear ok 
clear data_copy
Zthresh = 1.5; %threshold for rejecting artifacts
for i = 1:size(z_data,1)
    smoothed_abs_zdata(i,:) = conv(abs(z_data(i,:)),ones(300,1)/300,'same');
    ok(i,:) = smoothed_abs_zdata(i,:)<Zthresh; 
end
 data_copy = z_data; 
data_copy(~ok) = nan;

figure;
subplot(2,1,1);imagesc(z_data)
subplot(2,1,2);imagesc(data_copy)


%% re-reference data_clean 

%% re-referencing
cfg = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
data_eeg        = ft_preprocessing(cfg,data_eye);

% check whether this is doing what I think it is? 


%% 

data = data_eeg.trial{1}(1:61,:); 
data(~ok) = nan; 
data_eeg.trial{1}(1:61,:) = data;

%% redefine trials 
bl = 2; 
cfg = []; 
cfg = cfg_short_trials; 



if bl == 1
  
    cfg.trl = cfg.trl(1:end_of_block(bl),:); 
    
else 
    
    cfg.trl = cfg.trl(end_of_block(bl-1)+1:end_of_block(bl),:); 
    
end 

data_trial = ft_redefinetrial(cfg,data_eeg); 


%% save file 
file_path = '/Users/maria/desktop/data.continous_rdk/';
file_name = 'sub003_sess002_itis_intl.mat'; 
full = fullfile(file_path,file_name); 
save(full,'data_trial'); 

%% mark channels in which is a nan period in a trial to exclude them 




for i = 1:length(data_trial.trial)
    
    vec_good_channels(:,i) = sum(isnan(data_trial.trial{i}(:,:)),2) == 0;
    
    
end 


vec_good_channels(62:end,:) = 0; 
%% sort data 

cor_right{1} = find(data_trial.trialinfo == 24); 
cor_right{2} = find(data_trial.trialinfo == 25); 
cor_right{3} = find(data_trial.trialinfo == 26); 

cor_left{1} = find(data_trial.trialinfo == 34); 
cor_left{2} = find(data_trial.trialinfo == 35); 
cor_left{3} = find(data_trial.trialinfo == 36); 


%%  ft_timelock 


% highest coherence 


% right 
for i = 1:length(cor_right{1})
    tr = cor_right{1}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40); 
    data_cor_right_high(i,:) = data_trial.trial{tr}(40,:); 
   
end


% left 
for i = 1:length(cor_left{1})
    tr = cor_left{1}(i); 
    channels = zeros(65,1);
     channels(40) = channels(40) ; 
    data_cor_left_high(i,:) = data_trial.trial{tr}(40,:); 
   
end


mean_all_high = nanmean([data_cor_right_high;data_cor_left_high]); 

figure
subplot(3,1,1)
plot(mean_all_high,'b')
ylabel('mean voltage')
xlabel('msec')
title('high coherences')

%% middle coherences 


% right 
for i = 1:length(cor_right{2})
    tr = cor_right{2}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40) + vec_good_channels(40:tr); 
    data_cor_right_mid(i,:) = data_trial.trial{tr}(logical(channels),:); 
   
end


% left 
for i = 1:length(cor_left{2})
    tr = cor_left{2}(i); 
    channels = zeros(65,1);
     
    data_cor_left_mid(i,:) = data_trial.trial{tr}(40,:); 
   
end


mean_all_mid = nanmean([data_cor_right_mid;data_cor_left_mid]); 

subplot(3,1,2)
plot(mean_all_mid,'r')
ylabel('mean voltage')
xlabel('msec')
title('mid coherences')

%% low coherences 


% right 
for i = 1:length(cor_right{3})
    tr = cor_right{3}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40) + vec_good_channels(40,tr); 
    data_cor_right_low(i,:) = data_trial.trial{tr}(logical(channels),:); 
   
end


% left 
for i = 1:length(cor_left{3})
    tr = cor_left{3}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40) + vec_good_channels(40,tr); 
    data_cor_left_low(i,:) = data_trial.trial{tr}(logical(channels),:); 
   
end


mean_all_low = nanmean([data_cor_right_low;data_cor_left_low]); 

subplot(3,1,3)
plot(mean_all_low,'g')

ylabel('mean voltage')
xlabel('msec')
title('low coherences')

%% look at single channels across trials 




%% 
clear data_cor_right_high
clear data_cor_left_high 
clear data_cor_left_mid
clear data_cor_right_mid
clear data_cor_right_low
clear data_cor_left_low
