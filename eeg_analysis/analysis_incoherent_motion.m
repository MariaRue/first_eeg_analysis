% correlation of coherence values with EEG signal 
addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data 
addpath('/Users/maria/MATLAB-Drive/eeglab14_1_2b'); % to run eeglab to read in files and convert them 
addpath(genpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/converted_EEG_data')); % data location
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub001/behaviour');
ft_defaults % start fieldtrip 
eeglab % start eeglab
%% load in data 
clear all; 
close all; 
ft_defaults

 cfg                            = [];
 cfg.dataset                =  'sub001_sess017_eeg.set'; % your filename with file extension;
  cfg.trialdef.eventtype  = 'trigger'; 

 cfg.trialdef.eventvalue = [11]; % your event values
 cfg.trialdef.prestim    = 0;  % before stimulation (sec), only use positive num

 cfg.trialdef.poststim   = 300; % after stimulation (sec) , only use positive num

cfg                            = ft_definetrial(cfg);

data   = ft_preprocessing(cfg);
% downsample to 100hz 
cfg = []; 
cfg. resamplefs = 100; 
cfg. detrend = 'no'; 
cfg. demean = 'no'; 

data = ft_resampledata(cfg, data); 

%% %% re-referencing 
cfg = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
data_eeg        = ft_preprocessing(cfg,data);

%% filter data 

% bandpass filter the data 
cfg = []; 
% cfg.bpfilter = 'yes'; 
% cfg.bpfreq = [0.1, 30]; % problem with hp filter remains - not able to filter at 0.1hz what we want to do 
cfg.demean = 'yes';
cfg.detrend = 'yes';
% cfg.bpfilttype = 'firws';
%cfg.bpfiltord = 19;
%cfg.bpinstabilityfix = 'split'

cfg.hpfilter = 'yes';
cfg.hpfreq = 0.1; 
cfg.lpfilter = 'yes';
cfg.lpfreq = 30; % we played around with filtering to get rid of the weird artefacts we found, usually lowpass filter with 0.1 and ord 3 only 
cfg.lpfiltord = 3;
data_filtered = ft_preprocessing(cfg,data_eeg);
%% view data 
cfg = []; 
%cfg.viewmode = 'vertical'; 
cfg = ft_databrowser(cfg, data_eeg); 




%% get trigger info 

cfg = []; 
cfg.dataset = 'sub001_sess017_eeg.set';
cfg.trialdef.eventtype = 'trigger'; 
dummy = ft_definetrial(cfg); 

dummy.event


%% load behav data 
behav = load('sub001_sess017_behav.mat');

%% multiple regression 

% loop through each stimulus 

% sess19 stim counter 
% stim_count = 2; 
for i = 1:4
% for i = 1:4; 
   
stimulus = behav.S.coherence_frame{i}(1:30000); 
eeg_stream = data_filtered.trial{i};

% stim_count = stim_count + 1; 

% build EV with different time lags up to 1000ms out of stim stream and
% regress all of these with the EEG stream 
EVs{i} = zeros(size(eeg_stream,2),51); 

idx_end_stim = 0;

for l = 1:50
    
    EVs{i}(l:end,l) = stimulus(1:end-idx_end_stim); 
    
    idx_end_stim = idx_end_stim+1;
end

EVs{i}(:,51) = data_filtered.trial{i}(63,:);
% glmfit 
 
DV{i} = eeg_stream';
% DV{i} = DV(1001:end-1000,:);
% EVs{i} = EVs(1001 : end - 1000, :);

end 

%% betas for each block 

for i = 1:4
    X = EVs{i};
    pX = inv(X'*X)*X'; %pseudo inverse of EVs needed to calculate betas -
    % dont use pinv - takes too long - this is faster
    Y = DV{i}; % data
    beta{i} = pX*Y; % pseudo inverse of Evs times data gives betas
    
    
    % tstats
    
    tc = eye(size(X,2)); % we have 1000 different Evs here (for each lag one)
    tc(end+1,:) = sum(tc(:,1:end-1),2);
    
    prevar=diag(tc*pX*pX'*tc');
    R=eye(size(X,1))-X*pX; % this line of code does not work - would return a
    % 600 000 x 600 000 matrix - which is too big for mac memory!
    tR=trace(R);
    
    pe=pX*Y;
    cope=tc*pe;
    
    res=Y-X*pe;
    sigsq=sum(res.*res/tR);
    varcope=prevar*sigsq;
    
    tstat{i}=cope./sqrt(varcope);
    keyboard;
end

%% topoplot of betas 


data_beta = data_filtered; 

% betas in trial format for fieldtrip
    
     for i = 1:4
    
    data_beta.trial{i} = beta{i}(1:50,:)';
    
    data_beta.time{i} = data_filtered.time{i}(1:50);
     end 
     
     %%
     
b = beta{1}; % betas for first block 
Y = DV{1}; % original eeg data for first block 
E = EVs{1}; % Evs for first block
yey = E(:,51) * b(51,:); % predicted eye blinks 
ye = E(:,1:50) * b(1:50,:); % predicted data without eyeblinks 
Ywe = Y-yey; % original eeg - predicted eyeblinks 

% now compare Ywe with ye, 

figure 
plot(Ywe(:,5)); 
hold on 
plot(ye(:,5)); 


     


%%

cfg = [];                            
% cfg.xlim = [0.3 0.5];                
% cfg.zlim = [0 6e-14];                
cfg.layout = 'quickcap64.mat';
% cfg.trials = 1;
% 
%  cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_beta); colorbar


%% tstats topo

data_tstat = data_filtered;





% betas in trial format for fieldtrip
    
 for i = 4
    
    data_tstat.trial{i} = tstat{i}(1:50,:)';
    data_tstat.trial{i} = tstat{i}(1:50,:)';
    data_tstat.time{i} = data_filtered.time{i}(1:50);

 end 

 %%
     mydata =[];
     mydata.label = data_filtered.label;
     mydata.fsample = data_filtered.fsample;
     mydata.elec = data_filtered.elec;
     mydata.dimord = 'chan_time';
     mydata.time = 1;
    
     mydata.tstat = tstat{1}(52,:)';
   
     mydata.cfg = [];
     
     cfg = [];
     cfg.parameter = 'tstat';
     cfg.trial = 1; 
     cfg.zlim = [-6 6];
     cfg.layout = 'quickcap64.mat';
     figure; ft_topoplotTFR(cfg,mydata);
%%
cfg = [];                            
% cfg.xlim = [0.3 0.5];                
 cfg.zlim = [-3e-14 3e-14];                
cfg.layout = 'quickcap64.mat';
cfg.trials = 1;


  %cfg.parameter = 'none'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_tstat); colorbar

%% glm fit 

% build EVs 

% for i = 1:4; 

block = 1; 
   
stimulus = behav.S.coherence_frame{block}(1:30000); 
eeg_stream = data_filtered.trial{block};

% stim_count = stim_count + 1; 

% build EV with different time lags up to 1000ms out of stim stream and
% regress all of these with the EEG stream 
EVs = zeros(size(eeg_stream,2),50); 

idx_end_stim = 0;

for l = 1:50
    
    EVs(l:end,l) = stimulus(1:end-idx_end_stim); 
    
    idx_end_stim = idx_end_stim+1;
end

% glmfit 
 
DV = eeg_stream';
% DV{i} = DV(1001:end-1000,:);
% EVs{i} = EVs(1001 : end - 1000, :);

%EVs = EVs(:,10); 
%%

for i = 1:60
    [B{i},~,st{i}] = glmfit(EVs, DV(:,i),'normal','constant','off'); 
end

%% tstats

tstats = zeros(50,60);
for i = 1:50
    tstats(:,i) = st{i}.t;
end

%% 

for ch = 1 : 60
[C{ch},lag] = xcorr(data_filtered.trial{1}(ch),behav.S.coherence_frame{1},100);

end

%% autocorrelation stimulus 

[AC,lag] = xcorr(behav.S.coherence_frame{1},behav.S.coherence_frame{1},100);
