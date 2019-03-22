
subID = 9; % specify subID (is only number in file name)
save_name = sprintf('betas_sub%03.0f.mat',subID);


load(fullfile(EEGdatadir,save_name))% you need to specify EEGdatadir - path to the files you want to read in 


for r = 1:length(betas)-2
    figure;
    if tf_analysis
            plotmse(squeeze(betas{frequency}{r}(channel_ind,:,:)),2,time_idx(r).timebins);
              title(sprintf('Channel: %s Frequency: %dHz' ,chanlabel, freqlabel));
    else 
    plotmse(squeeze(betas{r}(channel_ind,:,:)),2,time_idx(r).timebins);
  %plot(time_idx(r).timebins,squeeze(betas{r}(channel_ind,:,:)));
    title(sprintf('Channel: %s' ,chanlabel));
    end
  xlabel(sprintf('Influence of %s on EEG at time (t+X) ms',time_idx(r).name));
    tidyfig;
end