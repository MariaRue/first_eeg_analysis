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
% load('data_03.mat');

%% load data JUST TO FIND OUT WHERE THE EPOCH BEGINS AND ENDS (in terms of samples)

cfg                        = [];
cfg.dataset                = 'sub003_sess001_fil001.set';
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

cfg = [];
cfg.dataset                 = 'sub003_sess001_fil001.set';
% cfg.trialdef.eventtype      = 'trigger';
% first and last samples

first_end_of_block = find(cfg_short_trials.trl(:,4)==210,1,'first');

cfg.trl                     = [cfg_short_trials.trl(1) cfg_short_trials.trl(first_end_of_block,2) 0];

% 
% filter 
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;

data_one_big_epoch          = ft_preprocessing(cfg);
%% downsample 
cfg = []; 
cfg.resamplefs = 150;
data_down = ft_resampledata(cfg,data_one_big_epoch);
%% regress out EOG channel 

X(:,1) = data_down.trial{1}(63,:); %Vertical EOG
X(:,1) = X(:,1) - mean(X(:,1)); %demean
X(:,2) = ones(size(X(:,1))); %constant term
for i = 1:61
    
    Y(i,:) = data_down.trial{1}(i,:)'; 

    betas(i,:) = glmfit(X,Y(i,:),'normal','constant','off');
    
end 

predYblink = betas(:,1)*X(:,1)';
imagesc(Y - predYblink);



%% find artifacts 

% find jumps 


cfg                    = [];

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'EEG';
cfg.artfctdef.zvalue.cutoff     = [8];
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0.05;
cfg.artfctdef.zvalue.fltpadding = 1;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfgeye, artifact_jump] = ft_artifact_zvalue(cfg,data_down);

% with this we are missing downward jumps... 


%%  detect muscles artifacts - this obviously only works BEFORE I lp filter the data - is that correct way of doing this? First find muscle artifact then filter? 
cfg = []; 
cfg.artfctdef.zvalue.channel = 'EEG';
  cfg.artfctdef.zvalue.cutoff      = 4;
  cfg.artfctdef.zvalue.trlpadding  = 0;
  cfg.artfctdef.zvalue.fltpadding  = 0;
  cfg.artfctdef.zvalue.artpadding  = 0.1;

  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter    = 'yes';
  cfg.artfctdef.zvalue.bpfreq      = [110 140]; % error message - butter: critical frequencies must be in (0 1) - why is this not working? 
  cfg.artfctdef.zvalue.bpfiltord   = 9;
  cfg.artfctdef.zvalue.bpfilttype  = 'but';
  cfg.artfctdef.zvalue.hilbert     = 'yes';
  cfg.artfctdef.zvalue.boxcar      = 0.2;

  % make the process interactive
  cfg.artfctdef.zvalue.interactive = 'yes';

  [cfgmuscle, artifact_muscle] = ft_artifact_zvalue(cfg,data_down);
  
  %% find heart - necessary? probably only in some participants who just have a very strong heart signal for some reason? 
  

cfg = []; 
   cfg.artfctdef.zvalue.channel     = 'EEG';
   cfg.artfctdef.zvalue.cutoff      = 4;
   cfg.artfctdef.zvalue.trlpadding  = 0;
   cfg.artfctdef.zvalue.artpadding  = 0.1;
   cfg.artfctdef.zvalue.fltpadding  = 0;

   % algorithmic parameters
   cfg.artfctdef.zvalue.bpfilter   = 'yes';
   cfg.artfctdef.zvalue.bpfilttype = 'but';
   cfg.artfctdef.zvalue.bpfreq     = [1 15];
   cfg.artfctdef.zvalue.bpfiltord  = 4;
   cfg.artfctdef.zvalue.hilbert    = 'yes';

   % feedback
   cfg.artfctdef.zvalue.interactive = 'yes';

   [cfg, artifact_EOG] = ft_artifact_zvalue(cfg);
   
   %% reject eye artifacts and padd with NaNs 
   cfg = []; 
   cfg.artfctdef.reject = 'nan'; 
   cfg.channel = 'EEG';
  cfg.artfctdef.jump.artifact = artifact_jump;
   data_no_artifacts  = ft_rejectartifact(cfg,data_down);
    
    %% look at data again and reject other artifacts manually 
cfg = [];
cfg.viewmode = 'vertical';

cfg.channel = 'EEG'; 


cfg = ft_databrowser(cfg,data_no_artifacts);

%% reject visual to identify bad channels 7


   cfg.artfctdef.reject = 'nan'; 
   cfg.channel = 'EEG';

data_clean = ft_rejectartifact(cfg,data_no_artifacts);


%% 
%% do component analysis on data_one_big_epoch to identify artifacts 

% output is 'comp'
% sampleinfo of 'comp' just the first and last samples of the epoch, as above

cfg        = [];
cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
cfg.channel = 'EEG';

comp = ft_componentanalysis(cfg, data_down); % ICA

%% plot components  
figure
cfg = [];
 cfg.component = 32;     % specify the component(s) that should be plotted
cfg.layout    = 'EasycapM1.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';

ft_topoplotIC(cfg, comp)
%% 

cfg = [];
% cfg.layout = 'quickcap64.mat'; % specify the layout file that should be used for plotting
%cfg.viewmode= 'vertical';
cfg.viewmode = 'component';
% cfg.channel = [1 2 3 7 9 15 38];

ft_databrowser(cfg, comp)



%% re-referance 
cfg = []; 
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.implicitref = 'LM';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel     = {'LM', 'RM'}; % the average of these channels is used as the new reference



data          = ft_preprocessing(cfg,data_clean);




%% redefine ICA'd data into short trials

cfg = [];
cfg.trl = cfg_short_trials.trl;
ft_redefinetrial