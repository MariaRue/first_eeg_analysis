EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
load(fullfile(EEGdir,'preprocessed_EEG_dat',sprintf('sub%03.0f_betas_kernel_reg.mat',20)));
%%

for r = 1:6
    figure;

    plotmse(squeeze(betas{r}(channel_ind,:,:)),2,time_idx(r).timebins);
  %plot(time_idx(r).timebins,squeeze(betas{r}(channel_ind,:,:)));
    title(sprintf('Channel: %s' ,chanlbCPZ));

  xlabel(sprintf('Influence of %s on EEG at time (t+X) ms',time_idx(r).name));
    tidyfig;
end
%%
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults

clear bs 
clear mean_b
ft_struct.time = time_idx(6).timebins;
bs = betas{1}(:,:,:,:);

% take the average across sessions
mean_b = nanmean(bs,4);
mean_b = nanmean(mean_b,3);
ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:);

%% 
cfg = [];
% cfg.xlim = [0.3 0.5];  % time limit
cfg.zlim = [-2 1];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,ft_struct); colorbar


