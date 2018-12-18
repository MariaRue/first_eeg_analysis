% this is a preprocessing for analysing data for the p300 signal 

addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data 
addpath('/Users/maria/MATLAB-Drive/eeglab14_1_2b'); % to run eeglab to read in files and convert them 
addpath(genpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/converted_EEG_data')); % data location
ft_defaults % start fieldtrip 
eeglab % start eeglab
%% first step is converting the curry files into fieldtrip compatible .set files
% via eeglab 

% pre-processsing pipeline for P300, not suitable for convolution - highpass filter
% at 0.1 removed, should be done for convolution because of skin potential
% that causes downward trend in eeg data, however for P300 analysis this is
% not so important because we cut out really small parts of the data 

% the cleaned data from this code was save in .mat files in the same folder
% and then analysed with the code p300_analysis.m 

%% try to read in events
% read in list with file names - list generated in terminal ls > ....txt,
% make sure only set files are in that list
fid = fopen('low_freq_sess.txt'); % list of files I want to append and analyse 
f = textscan(fid,'%s'); 

% loop through files extract coherent motion events and append files to
% each other 
num_sess = size(f{:},1); 

gap = 10000;

for i = 1 : num_sess

 cfg                            = [];
 cfg.dataset                =  f{1}{i}; % your filename with file extension;
 cfg.trialdef.eventtype  = 'trigger'; 

 cfg.trialdef.eventvalue = [25 35 50 75  125   135   150   175]; % your event values
 cfg.trialdef.prestim    = 0.5;  % before stimulation (sec), only use positive num

 cfg.trialdef.poststim   = 1.5; % after stimulation (sec) , only use positive num

cfg                            = ft_definetrial(cfg);

data{i}   = ft_preprocessing(cfg);

end

% fieldtrip has to be hacked otherwise it thinks that trials overlap
% because the sampleinfo is overlapping between some trials 
for i = 1 : num_sess - 1
    maxsamp = max(data{i}.sampleinfo(:)) + gap;
    data{i+1}.sampleinfo = data{i+1}.sampleinfo + maxsamp;
end


% append all data 
cfg = []; 
data_append = ft_appenddata(cfg,data{:}); 


%% re-referencing 
cfg = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
 data_eeg        = ft_preprocessing(cfg,data_append);
 %
 %% lowpass filter the data 
 cfg = []; 
 
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.lpfilter = 'yes';
cfg.lpfreq = 30; % we played around with filtering to get rid of the 
% weird artefacts we found, usually also highpass filter with 0.1 and ord 3  
cfg.lpfiltord = 3;
data_filtered = ft_preprocessing(cfg,data_eeg);


%% run ICA 
 
 % perform the independent component analysis (i.e., decompose the data)
cfg        = [];
cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB

comp = ft_componentanalysis(cfg, data_filtered); % ICA 
%%
% now plot the different components to decide which ones to remove from the
% data - we only removed eye movement components and components that looked
% clearly weird - started with topoplot because easier to identify
% eye movement components and outlayer components 
figure
cfg = [];
cfg.component = 1:24;       % specify the component(s) that should be plotted
cfg.layout    = 'quickcap64.mat'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';

ft_topoplotIC(cfg, comp)
%% browse through components to validate 

cfg = [];
cfg.layout = 'quickcap64.mat'; % specify the layout file that should be used for plotting
cfg.viewmode= 'vertical';
cfg.viewmode = 'component';
% cfg.channel = [1 2 3 7 9 15 38];
figure
ft_databrowser(cfg, comp)


 %% now reject components 
cfg = []; 
cfg.component = [7 12]; 

[data_final] = ft_rejectcomponent(cfg, comp, data_filtered);


%% look at data again 


cfg = [];
cfg.viewmode = 'vertical';


cfg = ft_databrowser(cfg,data_final);

%% below are additional methods to reject parts of the data that can be used to further clean the data
%% check for artifacts

cfg = []; 
cfg.channel = 'EEG';
cfg.method = 'summary';
[data_cleaned_more] = ft_rejectvisual(cfg,data_filtered); 
%% check for artifacts 
cfg = []; 
cfg.method = 'channel';
[data_clean] = ft_rejectvisual(cfg,epoched_data_filtered); 

%% check for artifacts 2
cfg = [];
cfg.viewmode = 'vertical';

cfg = ft_databrowser(cfg,data_eeg);

%% this function is used to reject artifacts 

data_cleaned = ft_rejectartifact(cfg,data_filtered);

%% check for artifacts (how Laurence does it)

cfg = []; 
cfg.channel = 'EEG';
cfg.method = 'summary';
[data_cleaned_more] = ft_rejectvisual(cfg,data_filtered); 

%% Pre-lim 300 analysis
% now separate for different coherences, pool left and right coherences for each coherence level 

coherence = unique(data_final.trialinfo); 
coherence = coherence(1:4);


for i = 1 : 4 % sort for coherences 
   
    idx_coh = data_final.trialinfo == coherence(i) | data_final.trialinfo == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_coherence{i} = ft_selectdata(cfg,data_final); 
    cfg = [];
average_ERP{i} = ft_timelockanalysis(cfg,data_coherence{i});

end

cfg = [];
%cfg.channel = {'P7';'P5';'P1';'PZ';'P2';'P3';'P4';'P6';'P8';'P07';'P03';'P0Z';'P04';'P08'};
cfg.baseline = [-0.4 -0.1];
cfg.baselinetype = 'absolute';
cfg.layout = 'quickcap64.mat';

figure; % i think this is wrong and I plotted these differently in the actual analysis 
hold on 
for i = 1:4
plot(average_ERP{1}.time,mean(average_ERP{i}.avg));
end
hold off
% 
% ft_singleplotER(cfg,average_ERP{3});


%% first attempts of data cleaning and analysis - not necessary 
% 
% % cfg = [];
% evs = EEG.event;
% toi = [25 35 50 75  125   135   150   175]; %triggers of interest
% pre_stim = 500; 
% post_stim = 1500;
% 
% events_of_interest = find(ismember([evs(:).type],toi));
% 
% count = 0;
% for i = events_of_interest
%     count = count+1;
%     trl(count,1) = evs(i).latency-pre_stim;
%     trl(count,2) = evs(i).latency+post_stim;
%     trl(count,3) = -pre_stim; 
% end
% 
% %now epoch the data
% % cfg.trl = trl;
% % epoched_data_filtered = ft_redefinetrial(cfg,data_filtered);
% 
% % steps above not needed because done at very start of script 
% 
% % now compute average response
% cfg = [];
% average_ERP = ft_timelockanalysis(cfg,data_filtered);
% 
% 
% % %% all data not event locked now same as above, because we had to event
% % lock data right at the start to get events into fieldtrip structure 
% % timelocked_data = ft_timelockanalysis(cfg, data_filtered);



% %% ICA with fieldtrip using fastica - suggested by Tom for identifying artifacts
% 
% % perform the independent component analysis (i.e., decompose the data)
% cfg        = [];
% cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
% 
% comp = ft_componentanalysis(cfg, data_cleaned_more); % ICA 
% %%
% % now plot the different components to decide which ones to remove from the
% % data 
% figure
% cfg = [];
% cfg.component = 1:20;       % specify the component(s) that should be plotted
% cfg.layout    = 'quickcap64.mat'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% 
% ft_topoplotIC(cfg, comp)
% %%
% 
% cfg = [];
% cfg.layout = 'quickcap64.mat'; % specify the layout file that should be used for plotting
% cfg.viewmode='vertical'
% cfg.viewmode = 'component';
% figure
% ft_databrowser(cfg, comp)
% 
% 
% 
% 
% %% 
% 
% %CP4, F4, POZ, FCZ, T8, P2, (PC5 or PC3), P3  
% 
% % run ft_select and remove bad channels then run ICA again 
% remove= [{'P4'}; {'F4'}; {'POZ'}; {'FCZ'}; {'T8'}; {'P2'}; {'PO5'}; {'PO3'}; {'P3'}];
% 
% 
% cfg = [];
% cfg.channel = data_cleaned_more.label(~ismember(data_cleaned_more.label, remove))
% data_channels_removed = ft_selectdata(cfg,data_cleaned_more);
% 
% 
% %% repeat ICA
% 
% 
% % perform the independent component analysis (i.e., decompose the data)
% cfg        = [];
% cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
% 
% comp = ft_componentanalysis(cfg, data_channels_removed); % ICA 
% %% 
% cfg = [];
% cfg.layout = 'quickcap64.mat'; % specify the layout file that should be used for plotting
% cfg.viewmode='vertical'
% cfg.viewmode = 'component';
% figure
% ft_databrowser(cfg, comp)
% 
% %%
% % now plot the different components to decide which ones to remove from the
% % data 
% figure
% cfg = [];
% cfg.component = 1:20;       % specify the component(s) that should be plotted
% cfg.layout    = 'quickcap64.mat'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% 
% ft_topoplotIC(cfg, comp)
% 
% %% now reject components 
% cfg = []; 
% cfg.component = [3 5 7 15 43]; 
% 
% [data_final] = ft_rejectcomponent(cfg, comp, data_cleaned_more);
% 
% 
% %% check final data 
% cfg = [];
% cfg.viewmode = 'vertical';
% 
% cfg = ft_databrowser(cfg,data_final);
% 
% %% average over parietal electrodes 
% cfg = [];
% cfg.channel = {'P7';'P5';'P1';'PZ';'P2';'P3';'P4';'P6';'P8';'P07';'P03';'P0Z';'P04';'P08'};
% cfg.baseline = [-0.4 -0.1];
% cfg.baselinetype = 'absolute';
% cfg.layout = 'quickcap64.mat';
% figure; ft_singleplotER(cfg,data_final);
% % use interactive mode to plot scalp topography.



