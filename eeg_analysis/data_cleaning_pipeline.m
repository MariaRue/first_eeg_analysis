% This script is for cleaning up EEG data from the continuous random dot
% experiment. Raw data recorded with curry have to be loaded into eeglab
% under file - import data - using EEGlab functions and plugins - load
% curry data (that is an additional package that needs to be downloaded to
% load in curry data) Data has to be saved in eeglab format which can then
% be uploaded to fieldtrip

%% only for transforming raw data into format that can be read by fieldtrip
addpath('/Users/maria/MATLAB-Drive/eeglab14_1_2b/')
eeglab

%% start fieldtrip
addpath('/Users/maria/MATLAB-Drive/fieldtrip-master/')
ft_defaults;

%% path to data
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub003/EEG/');

%%
load('data_03.mat');
%% load data

cfg                            = [];
cfg.dataset                = 'sub003_sess001_fil001.set';
cfg.trialdef.eventtype  = 'trigger';

cfg.trialdef.eventvalue = [24]; % your event values
cfg.trialdef.prestim    = 1;  % before stimulation (sec), only use positive num

cfg.trialdef.poststim   = 3; % after stimulation (sec) , only use positive num

cfg                            = ft_definetrial(cfg);

data   = ft_preprocessing(cfg);

%% how to take the entire epoch as one trial

cfg = [];
cfg.dataset                 = 'sub003_sess001_fil001.set';
% cfg.trialdef.eventtype      = 'trigger';
cfg.trl                     = [data.sampleinfo([1 end]) 0];

% cfg                         = ft_definetrial(cfg);
data_one_big_epoch          = ft_preprocessing(cfg);

%% re-referencing
cfg = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference
data_eeg        = ft_preprocessing(cfg,data);
%% downsample
cfg = [];
cfg.resamplefs = 100;
data_down = ft_resampledata(cfg,data_eeg);
%% take a look at individual trials

for i = 1:20
    plot(data_down.trial{i}(40,:));
    
    pause
    
end

for i = 1:308
    cpz(i,:) = data_down.trial{i}(40,:);
end

figure
imagesc(cpz)
colorbar

figure
plot(mean(cpz))

%%
for i = 1:308
    eye_blinks(i,:) = data_filtered.trial{i}(63,:);
end

figure
imagesc(eye_blinks)
colorbar

%%
cfg = [];
cfg.viewmode = 'vertical';

%cfg.ylim = 'maxmin';


cfg = ft_databrowser(cfg,data_down);

%% this function is used to reject artifacts

data_cleaned = ft_rejectartifact(cfg,data_filtered);



%%
cfg = [];

%cfg.demean = 'yes';
cfg.detrend = 'yes';
%cfg.lpfilter = 'yes';
cfg.lpfreq = 40; % we played around with filtering to get rid of the
% weird artefacts we found, usually also highpass filter with 0.1 and ord 3
% cfg.lpfiltord = 3;
data_filtered = ft_preprocessing(cfg,data_down);
%% check for artifacts
cfg = [];
cfg.channel = 'EEG';
[data_clean] = ft_rejectvisual(cfg,data_filtered);

%% run ICA

% perform the independent component analysis (i.e., decompose the data)
cfg        = [];
cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
cfg.channel = 'EEG';

comp = ft_componentanalysis(cfg, data_filtered); % ICA

%%
% now plot the different components to decide which ones to remove from the
% data - we only removed eye movement components and components that looked
% clearly weird - started with topoplot because easier to identify
% eye movement components and outlayer components
figure
cfg = [];
cfg.component = [37 48 52 59 61 63];       % specify the component(s) that should be plotted
cfg.layout    = 'quickcap64.mat'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';

ft_topoplotIC(cfg, comp)

%% browse through components to validate

cfg = [];
% cfg.layout = 'quickcap64.mat'; % specify the layout file that should be used for plotting
%cfg.viewmode= 'vertical';
cfg.viewmode = 'component';
% cfg.channel = [1 2 3 7 9 15 38];

ft_databrowser(cfg, comp)


%%

for i = [37, 48, 52, 59, 61, 63]
    for l = 1:308
        comp_vec(l,:) = comp.trial{l}(i,:);
    end
    figure
    subplot(3,1,1:2)
    imagesc(comp_vec);
    title(i)
    colorbar
    subplot(3,1,3)
    plot(mean(comp_vec));
    
    
end

