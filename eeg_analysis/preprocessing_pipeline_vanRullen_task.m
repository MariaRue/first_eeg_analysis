% this is a script to preprocess eeg data recorded from LH (subject 1) to
% check whether we can reproduce the results from vanRullen and MacDonald
% (2012 - 'Perceptual Echoes at 10Hz in the human brain'). 
%% add eeglab and fieldtrip paths and start both tbs 
addpath('/Users/maria/MATLAB-Drive/fieldtrip'); % fieldtrip tool box to analyse data 
addpath('/Users/maria/MATLAB-Drive/eeglab14_1_2b'); % to run eeglab to read in files and convert them 

ft_defaults % start fieldtrip 
eeglab % start eeglab

% path to data 

% behaviour 
addpath(genpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub001/behaviour')); % data location

% eeg (where transformed eeg data is saved, curry eeg data has to be transformed with eeglab into a set file, this can be read by fieldtrip) 
addpath(genpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/converted_EEG_data')); % data location
%% put eeglab data in fieldtrip structure sorted for trials. Each run of the 
% perceptual_echoes function from Laurence produces 45 trials, trigger for
% start of a trial = 4, each trial is 10secs long 

cfg = []; 
cfg.dataset = 'sub01_perceptual_echoes_2.set'; 
cfg.trialdef.eventtype = 'trigger'; 
cfg.trialdef.eventvalue = 4; 
cfg.trialdef.prestim = 0;
cfg.trialdef.poststim = 10; 
cfg = ft_definetrial(cfg); 

data = ft_preprocessing(cfg); 

%% re-referencing 
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
%cfg.bpinstabilityfix = 'split';
cfg.lpfilter = 'yes';
cfg.lpfreq = 50; % we played around with filtering to get rid of the weird artefacts we found, usually lowpass filter with 0.1 and ord 3 only 
cfg.lpfiltord = 3;
cfg.hpfilter = 'yes';
cfg.hpfreq = 1; 
cfg.hpfiltord = 3; 
data_filtered = ft_preprocessing(cfg,data_eeg);
%% view data 
cfg = []; 
%cfg.viewmode = 'vertical'; 
cfg = ft_databrowser(cfg, data_filtered); 

%% read in behavioural data 
percept_behav1 = load('LH03072018_test2.mat');
%% we only have 44 trials in eeg data, but there are 45 trials in behaviour.
% I assume that the first trial is missing but I am not sure and haven't 
% found a way yet to determine which trial is missing in the EEG data - 
% this is true for both sessions we recorded for this task - VEP only
% visible if 2nd trial of stimulus data is matched up with 1st trial of EEG
% data!!!! 
%% calculate cross correlation between stimulus and EEG for each channel for 1sec lag in both directions 
count_stim  = 2; 
C = zeros(65,2001,44);
for t = 1:44
for i = 1:65
[C(i,:,t),lags{i}] = xcorr(data_filtered.trial{t}(i,:),percept_behav1.trialVariables(count_stim).stimstream,'coeff',1000); % normalisation necessary? 
end
count_stim = count_stim+1; 
end

%% take mean average for channes over occipital/parietal c 
mean_IRF = zeros(11,2001);
count = 1; 
for ch = 50:61
    
    mean_IRF(count,:) = mean(C(ch,:,:),3);
    count = count + 1; 
end 


%% mean across channels 
figure;
mean_channels = mean(mean_IRF,1);

plot(lags{1},mean_channels);
%% plot results 
for i = 1:11
subplot(3,4,i)
plot(lags{1},mean_IRF(i,:),'b');
hold on;
end
%% 
 figure; 
for tr = 1:44
    
   
    hold on
    plot(C(55,:,tr)); 
    
    
end

%% plot correlation topography 

% transfer correlation data back into a fieldtrip structure 

data_corr = data_filtered; 

% correlation in trial format for fieldtrip 

for tr = 1:44
    
    data_corr.trial{tr} = C(:,:,tr);
    data_corr.time{tr} = data_filtered.time{1}(1:2001);
 
end

cfg = [];                            
% cfg.xlim = [0.3 0.5];                
% cfg.zlim = [0 6e-14];                
cfg.layout = 'quickcap64.mat';          
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_corr); colorbar
%% check whether GLM gives same results as calculated for correlation on single trial basis 


% build EV with different time lags up to 1000ms out of stim stream for
% each trial 
for tr = 1:44
EVs{tr} = zeros(size(data.time{1},2),1000); 
end
count_stim = 2; %(first eeg trial has to bee regressed with second stimulus, because first trial was not recorde)

for tr = 1:44
    idx_end_stim = 0;
for i = 1:1001
    if i == 1
        
       EVs{tr}(i:end,i) = ones(1,size(data.time{1},2));   
        
    end
    
    EVs{tr}(i:end,i) = percept_behav1.trialVariables(count_stim).stimstream(1:end-idx_end_stim); 
    
    idx_end_stim = idx_end_stim+1;
    
end
count_stim = count_stim+1;
end
%% fit glm for each trial 
beta_per_trial = zeros(1001,65,44); 
for tr = 1:44
X = EVs{tr};
pX = inv(X'*X)*X'; %pseudo inverse of EVs needed to calculate betas - 
% dont use pinv - takes too long - this is faster 
Y = data_filtered.trial{tr}'; % data 
beta_per_trial(:,:,tr) = pX*Y; % pseudo inverse of Evs times data gives betas 
end

mean_beta = mean(beta_per_trial,3); 
%% topoplot o_beta_per_trial - get beta_per_trial back in fieldtrip structure 


% transfer beta data back into a fieldtrip structure 

data_beta = data_filtered; 

% betas in trial format for fieldtrip 

for tr = 1:44
    
    data_beta.trial{tr} = beta_per_trial(:,:,tr)';
    data_beta.time{tr} = data_filtered.time{1}(1:1001);
 
end

cfg = [];                            
% cfg.xlim = [0.3 0.5];                
% cfg.zlim = [0 6e-14];                
cfg.layout = 'quickcap64.mat';          
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_beta); colorbar
%% now calculate the tstats for this 


% define a contrast matrix (in this case it is an identity matrix from number of EVs
% (but doesn't have to be one!!! In this case it is one because we want to
% test whether all betas from each EV are significantly different to 0 

tc = eye(size(EVs{1},2)); % we have 1000 different Evs here (for each lag one)
tstat = zeros(size(EVs{1},2),65,44);
for tr = 1:44
    
    
 Y = data_filtered.trial{tr}';
 X = EVs{tr};
 pX = inv(X'*X)*X';
prevar=diag(tc*pX*pX'*tc');

id = eye(size(X,1)); 
trm = X*pX; 
R=id - trm; % this line of code does not work - would return a 
% 600 000 x 600 000 matrix - which is too big for mac memory!
tR=trace(R);

pe=pX*Y;
cope=tc*pe;

res=Y-X*pe;
sigsq=sum(res.*res/tR);
varcope=prevar*sigsq;


tstat(:,:,tr)=cope./sqrt(varcope);
end

%% topoplot of tstat- get tstats per trial back in fieldtrip structure 


% transfer beta data back into a fieldtrip structure 

data_tstat = data_filtered; 

% betas in trial format for fieldtrip 

for tr = 1:44
    
    data_tstat.trial{tr} = tstat(:,:,tr)';
    data_tstat.time{tr} = data_filtered.time{1}(1:1001);
 
end

cfg = [];                            
% cfg.xlim = [0.3 0.5];                
% cfg.zlim = [0 6e-14];                
cfg.layout = 'quickcap64.mat';          
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_tstat); colorbar



%% simulate data from random noise input 
% functions taken from VanRullen CCN talk 2017
% (https://www.youtube.com/watch?v=r1g7kj0HjZo) 

deltaT = 10; % processesing delay between layer 1 and layer 2 (see below)
tau = 15; 
%input = rand(100,1);

num_trial= 1000;

y = zeros(size(input)); 
x = zeros(size(input));

trial = zeros(num_trial,199);


passbandfreq    = 10;
stopbandfreq    = 15;
passrip         = 1;
stopbandatten   = 60;
framerate       = 60; 
coherence_sd    = 0.05;
ft = designfilt('lowpassfir', 'PassbandFrequency', passbandfreq, 'StopbandFrequency', stopbandfreq,...
    'PassbandRipple', passrip, 'StopbandAttenuation', stopbandatten, 'SampleRate', framerate);


for tr = 1:num_trial % loop through simulations 
input = randn(1000,1); % generate random noise 
input = filter(ft,input); 
y = zeros(size(input)); % values from function predicting next sample input (layer 1)
x = zeros(size(input)); % values from function calculating prediction error from y and sample input (layer 2)

for t =   1:length(input) % loop through input stream 
    
    
     if t - deltaT > 0 


    x(t) = input(t) - y(t-deltaT);
    
    y(t+1) = y(t) + x(t-deltaT)/tau-(y(t)/200);
    

     else 
    x(t) = input(t) - y(t);
    
    y(t+1) = y(t) + x(t)/tau;
% 
    end
%     
end
stimmean = (input - mean(input))/std(input);
ymean = (y - mean(y)) / std(y);

 if tr == 1
     
     figure 
     plot(stimmean)
     hold on 
     plot(ymean)
     hold off
     xlabel('time')
     title('example white noise trial')
 end
     

 
[C,lag] = xcorr(ymean(1:100), stimmean(1:100),'coeff');


trial(tr,:) = C;


end

meantr = mean(trial);

figure;
plot(lag,meantr); 
xlabel('lag');
ylabel('correlation');
title('white noise');


%% simulate data with task stim 

deltaT = 20; 
tau = 30; 
input = rand(1000,1);

num_trial= 45;

y = zeros(size(input)); 
x = zeros(size(input));

trial = zeros(num_trial,2001);

for tr = 1:num_trial
input = percept_behav1.trialVariables(tr).stimstream;

y = zeros(size(input)); 
x = zeros(size(input));
for t =   1:length(input)
    
    
     if t - deltaT > 0 


    x(t) = input(t) - y(t-deltaT);
    
    y(t+1) = y(t) + x(t-deltaT)/tau-(y(t)/200);
    

     else 
    x(t) = input(t) - y(t);
    
    y(t+1) = y(t) + x(t)/tau;
% 
    end
%     
end
stimmean = (input - mean(input))/std(input);
ymean = (y - mean(y)) / std(y);


 if tr == 1
     
     figure 
     plot(stimmean)
     hold on 
     plot(ymean)
     hold off
     xlabel('time')
     title('example stim trial')
 end
 
 
[C,lag] = xcorr(stimmean, ymean,1000);


trial(tr,:) = C;


end

meantr = mean(trial);

figure;
plot(lag,meantr); 
xlabel('lag');
ylabel('correlation');
title('stimulus');


%% try cross correlation across several trials combined as continuous data, choose an arbitrary amoutn of time between indiviual 'trials' and fill with zeros 


total_stim_samples = 45 .* 10000; 
time_bet_trials = 1000; 
total_time_bet_trials = 100 .* 44; 

total_stimstream = zeros(1,total_stim_samples+total_time_bet_trials); 


% fill total_stimstream with stimuli from perceptual echoes session 
start_idx = 1; 
end_idx = 10000;

for tr = 1:45
    
    total_stimstream(1,start_idx : end_idx) = percept_behav1.trialVariables(tr).stimstream; 
    start_idx = end_idx + time_bet_trials; 
    end_idx = start_idx + 10000-1; 
    
end


% simulate brain response for reach trial and then fill up with
deltaT = 20; 
tau = 30; 
Y = zeros(45,10000);
for tr = 1:45
input = percept_behav1.trialVariables(tr).stimstream;

y = zeros(size(input)); 
x = zeros(size(input));
for t =   1:length(input)
    
    
     if t - deltaT > 0 


    x(t) = input(t) - y(t-deltaT);
    
    y(t+1) = y(t) + x(t-deltaT)/tau-(y(t)/200);
    

     else 
    x(t) = input(t) - y(t);
    
    y(t+1) = y(t) + x(t)/tau;
% 
    end
%     
end

Y(tr,:) = y(1:end-1); 
end

% fill up vector with all simulated brain signals combined 

total_brain_sim = zeros(1,total_stim_samples+total_time_bet_trials); % has to be same length as total_stimstream

% fill total_brain_sim with simulated brain signals in Y 
start_idx = 1; 
end_idx = 10000;

for tr = 1:45
    
    total_brain_sim(1,start_idx : end_idx) = Y(tr,:); 
    start_idx = end_idx + time_bet_trials; 
    end_idx = start_idx + 10000-1; 
    
end

% demean stim stream and brain sim stream 

dmean_stim = total_stimstream - mean(total_stimstream); 
dmean_brain = total_brain_sim - mean(total_brain_sim); 

% correlate both 

[C_stream_sim, lag] = xcorr(dmean_brain,dmean_stim,'coeff'); 

figure; 
plot(lag,C_stream_sim)
xlabel('lag ms')
ylabel('correlation')
title('correlation of stream with all trials with stream of all simulated signals, not segemented for trials')
%% correlate entire eeg signal with all trials (not segmented for trials) 

 cfg                            = [];
 cfg.dataset                =  'sub01_perceptual_echoes_2.set'; % your filename with file extension;
%   cfg.trialdef.eventtype  = 'trigger'; 
% 
%  cfg.trialdef.eventvalue = [11:19]; % your event values
%  cfg.trialdef.prestim    = 0;  % before stimulation (sec), only use positive num
% 
%  cfg.trialdef.poststim   = 36; % after stimulation (sec) , only use positive num
% 
% cfg                            = ft_definetrial(cfg);

data_continious   = ft_preprocessing(cfg);

%% re-referencing 
cfg = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
data_continious       = ft_preprocessing(cfg,data_continious);

%% filter data 

% bandpass filter the data 
cfg = []; 
% cfg.bpfilter = 'yes'; 
% cfg.bpfreq = [0.1, 30]; % problem with hp filter remains - not able to filter at 0.1hz what we want to do 
cfg.demean = 'yes';
cfg.detrend = 'yes';
% cfg.bpfilttype = 'firws';
%cfg.bpfiltord = 19;
%cfg.bpinstabilityfix = 'split';
cfg.lpfilter = 'yes';
cfg.lpfreq = 50; % we played around with filtering to get rid of the weird artefacts we found, usually lowpass filter with 0.1 and ord 3 only 
cfg.lpfiltord = 3;
cfg.hpfilter = 'yes'; 
cfg.hpfreq = 0.1; 
cfg.hpfiltord = 3;
data_continious = ft_preprocessing(cfg,data_continious);
%%
cfg = []; 
%cfg.viewmode = 'vertical'; 
cfg = ft_databrowser(cfg, data_continious); 
%% get trigger info 

cfg = []; 
cfg.dataset = 'sub01_perceptual_echoes_2.set';
cfg.trialdef.eventtype = 'trigger'; 
dummy = ft_definetrial(cfg); 

dummy.event


% save event values and sample info in matrices (can be more easily
% accessed) 
for i = 1 : 5064
 
  event_val(i) = dummy.event(i).value; 
  event_sample(i) = dummy.event(i).sample;
    
end

%% restructre eeg signals and stim info so that the whole time course of the 
% EEG signal can be visualised. 

% get length of samples between end of trial/beginng of trial and beginning
% of trial/end of trial 
idx = event_val == 4 | event_val == 5; % 4 trigger for start of trial, 5 trigger for end of trial 
diff_sample = diff(event_sample(idx)); 
diff_sample = diff_sample(2:end); % first entry is end of trial trigger to
% first trial trigger of trial 2 (trial 1 has been missed because recording 
% button has been pressed too late)
%%
% recordings started after first trial has already started. We now need to
% find the first trial start trigger that belongs to trial 2 in the stim
% info structure 

[event_code,first_trial_start,~] = unique(event_val,'first'); 
first_trial_sample_idx = event_sample(first_trial_start(1)); 


% get end time of last trial 
[event_code,last_trial_end,~] = unique(event_val,'last'); 
last_trial_sample_idx = event_sample(last_trial_end(2)); 
%%
% EEG stream 
eeg_stream = data_continious.trial{1}(:,first_trial_sample_idx:last_trial_sample_idx);



%%
% now put togehter the stimstream 
stimstream = zeros(1,size(eeg_stream,2)); 

start_idx = 1; 
end_idx = 10000; 
diff_vec_count = 1; 
for tr = 2:45
    
    stimstream(start_idx:end_idx) = percept_behav1.trialVariables(tr).stimstream;
    if tr < 45
    % calculate how much longer eeg signal is than stim signal 
    diff_stim_eeg = diff_sample(diff_vec_count) - 10000; 
    
    diff_vec_count = diff_vec_count+1; 
    
    % now calculate total number of samples between stimuli - diff_stim_eeg
    % + difference between last sample of prev trial until new trial 
   
    between_stim = diff_stim_eeg + diff_sample(diff_vec_count); 
    diff_vec_count = diff_vec_count+1; 
    
    start_idx = end_idx + between_stim+1; 
    end_idx = start_idx + 10000-1; 
  
    end
end
%% make eegstream values 0 where stimstream is 0 
idx_0 = stimstream == 0; 
for ch = 1:65
eeg_stream(ch,idx_0) = 0; 
end

% now demean eeg_stream and stimstream 
dmean_eegstream = zeros(size(eeg_stream));
for ch = 1:65
dmean_eegstream(ch,:) = eeg_stream(ch,:) - mean(eeg_stream(ch,:));  
end
dmean_stimstream = stimstream - mean(stimstream); 

C_eeg_stim = zeros(65,2.*length(dmean_stimstream)-1); 
% C_eeg_stim = zeros(65,2001); 
for ch = 1:65 % correlation for each channel
[C_eeg_stim(ch,:), lag] = xcorr(dmean_eegstream(1,:),dmean_stimstream,'coeff');
end
%%  plot correlation 
figure;
count_plot = 1;
for ch = 50 : 60
   
    subplot(4,3,count_plot)
    plot(lag,C_eeg_stim(ch,:));
    count_plot = count_plot + 1; 
end
subplot(4,3,1)
title('correlation across entire experiment')
ylabel('correlation')
xlabel('lag ms')

%% autocorrelate stim stream with itself 

[AC, Alag] = xcorr(dmean_stimstream,'coeff'); 
% it seems like this shows only correlation between stim on/stim off,
% therefore we try multiple regression 

% single trial stim
[AC, Alag] = xcorr(percept_behav1.trialVariables(1).stimstream);

%% multiple regression 

% build EV with different time lags up to 1000ms out of stim stream and
% regress all of these with the EEG stream 
EVs = zeros(size(eeg_stream,2),1000); 

idx_end_stim = 0;

for i = 1:1000
    
    EVs(i:end,i) = dmean_stimstream(1:end-idx_end_stim); 
    
    idx_end_stim = idx_end_stim+1;
end

% glmfit 
DV = dmean_eegstream';
DV = DV(1001:end-1000,:);
EVs = EVs(1001 : end - 1000, :);

%% 

for i = 1:60
[B{i},~,st{i}] = glmfit(EVs, DV(:,i)); 

end

%% 
figure 

for i = 1:6
subplot(3,2,i)

plot(B{i})

end

for i = 1:6
    tt(:,i) = st{i}.t;
end

figure 

for i = 1:6
subplot(3,2,i)

plot(tt(:,i))

end

%%
addpath('/Users/maria/MATLAB-Drive/MATLAB/fun/');
%%
% try it with ols in folder fun
tc = eye(1000); 
[cope,varcope,tstat,F,Fdof]=ols(DV,EVs,tc);

%% 
X = EVs;
pX = inv(X'*X)*X'; %pseudo inverse of EVs needed to calculate betas - 
% dont use pinv - takes too long - this is faster 
Y = DV; % data 
beta = pX*Y; % pseudo inverse of Evs times data gives betas 


% now calculate tstats, this is the contrast matrix times the betas divided
% by std of this (code from ols function from TB and LH)

% define a contrast matrix (in this case it is an identity matrix from number of EVs
% (but doesn't have to be one!!! In this case it is one because we want to
% test whether all betas from each EV are significantly different to 0 

tc = eye(size(X,2)); % we have 1000 different Evs here (for each lag one)

prevar=diag(tc*pX*pX'*tc');
R=eye(size(X,1))-X*pX; % this line of code does not work - would return a 
% 600 000 x 600 000 matrix - which is too big for mac memory!
tR=trace(R);

pe=pX*Y;
cope=tc*pe;

res=Y-X*pe;
sigsq=sum(res.*res/tR);
varcope=prevar*sigsq;

tstat=cope./sqrt(varcope);

%% topoplot of tstat- get tstats per trial back in fieldtrip structure 


% transfer beta data back into a fieldtrip structure 

data_beta_all = data_filtered; 

% betas in trial format for fieldtrip 


    
    data_beta_all.trial{1} = beta;
    data_beta_all.time{1} = data_filtered.time{1}(1:1001);
 


cfg = [];                            
% cfg.xlim = [0.3 0.5];                
% cfg.zlim = [0 6e-14];                
cfg.layout = 'quickcap64.mat';          
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,data_beta_all); colorbar


