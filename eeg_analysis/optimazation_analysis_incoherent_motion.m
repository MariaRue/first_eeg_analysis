% This script is for trying gradient descent as developed in the
% simulations with an OU process to recover the integration kernel of the
% brain during the incoherent motion periods 

addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data 
addpath('/Users/maria/MATLAB-Drive/eeglab14_1_2b'); % to run eeglab to read in files and convert them 
addpath(genpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/converted_EEG_data')); % data location
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub001/behaviour');
addpath(genpath('/Users/maria/Matlab-Drive/Matlab/continous_rdk/continuous_rdk_simulations'));
ft_defaults % start fieldtrip 

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
cfg.viewmode = 'vertical'; 
cfg = ft_databrowser(cfg, filtered); 

%% load behav data 
behav = load('sub001_sess017_behav.mat');

%% 

stimulus = behav.S.coherence_frame{1}(1:29999); 
eeg_data = data_filtered.trial{1}; 

%  define bounds of parameters lambda, offset, amplitude

bounds_l = [-1 0]; % lambda
bounds_o = [0 900]; % offset
bounds_A = [0 10]; % amplitude

% generate 10 equally spaced variables within these bounds
lambda_grid = logspace(bounds_l(1),bounds_l(2),10);
offset_grid = round(linspace(bounds_o(1),bounds_o(2),10));
amplitude_grid = linspace(bounds_A(1),bounds_A(2),10);


sse_org = 10000000000000000000000000000000000000000; % initial value new calculated cost function is compared to

% loop through EEG channels 
for ch = 1:60
% do a grid search for each EEG channel and save that data in a big
% structure 

disp(ch);
 % loop through all possible parameter combinations
                for l_test = 1:length(lambda_grid)
                    
                    for off_test = 1:length(offset_grid)
                        
                        for a_test = 1:length(amplitude_grid)
                            
                            p0(1) = lambda_grid(l_test);
                            p0(2)  = offset_grid(off_test);
                            p0(3) = amplitude_grid(a_test);
                            
                            sse = oueval(p0,stimulus,eeg_data(ch,:),0.1);
                        
                            cost_channel(ch).sse(l_test,a_test,off_test) = sse;
%                         if p0(1) == lambda_grid(4) && p0(2) == offset_grid(3) && p0(3) == amplitude_grid(6)
%                             keyboard; 
%                         end
                         

                            if sse < sse_org % save best estimates of p from grid search
%                                 disp(p0)
%                                 disp(sse)

                                pstart{ch} = p0;
                                sse_org = sse;
                                
                               
                            end
                            
                         
                        end % loop through amplitudes
                        
                    end % loop through offsets
                    
                end % loop through lambdas

% now do the true parameter recovery with fminsearch 



end % loop through channels 


%% plot energy maps for each offset and each channel 


for ch = 1:60 
     figure (ch)
     sub_index = 0; 
    for off = 1:10
    sub_index = sub_index + 1; 
    
   
    subplot(5,2,sub_index) 
    imagesc(cost_channel(ch).sse(:,:,off));
    
    if sub_index == 1; 
        ylabel('amplitude value')
        xlabel('lambda value') 
        title(['channel ', num2str(ch), ' off= ', num2str(offset_grid(off))])
    else 
        title(['off= ', num2str(offset_grid(off))])
    end 
    end % loop throuch different offset grids 
end % loop through channels 

%% run fminserach for all channels 
for ch = 1:60 
    disp(ch);
    if size(pstart{ch},2) == 3
     fun = @(p)oueval(p,stimulus,eeg_data(ch,:),0.1); % this is the correct cost function that works
                [pnew] = fminsearch(fun,pstart{ch});
                keyboard;
                pnew(1) = -abs(pnew(1));
                estimates(ch,:) = pnew;
                % pguess{noise,off}(:,rep,l) = pstart; 
    else 
        estimates(ch,:) = zeros(1,3);
    end
    
end % loop through channels 
