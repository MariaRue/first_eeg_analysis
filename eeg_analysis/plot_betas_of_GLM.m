% EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
% load(fullfile(EEGdir,'preprocessed_EEG_dat',sprintf('sub%03.0f_betas_kernel_reg.mat',20)));

scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
%EEGdir = '/Volumes/LaCie/data/preprocessed_EEG_dat/';
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','preprocessed_EEG_dat');
addpath('/Users/maria/Documents/matlab/spm12');
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults

% path to hard drive /Volumes/LaCie
subj_list = [16,18:21,24,26,32];
%%
fig_Id = 1;
figure_counter = 0;
for i = 1:length(subj_list)
    subID = subj_list(i);
    load(fullfile(EEGdir, sprintf('sub%03.0f_betas_kernel_reg.mat',subID)))
    
    for r = 1:4
        subplot(4,4,fig_Id)
        plotmse(squeeze(betas{r}(channel_ind,:,:)),2,time_idx(r).timebins);
        %plot(time_idx(r).timebins,squeeze(betas{r}(channel_ind,:,:)));
        
        
        xlabel(sprintf('Influence of %s on EEG at time (t+X) ms',time_idx(r).name));
        tidyfig;
        if r ==1
            title(sprintf('Subject: %3.0f Channel: %s',subID, chanlbCPZ))
            
        else
            title(sprintf('Channel: %s' ,chanlbCPZ));
            
            
        end
        
        if fig_Id == 16
            fig_Id = 1;
            
            
            figure_counter = figure_counter + 1;
            if figure_counter == 1
                savefig(fullfile(EEGdir,'4_subjs_kernel_reg1.fig'))
                figure;
            else
                savefig(fullfile(EEGdir,'4_subjs_kernel_reg2.fig'))
                
            end
            
            
        else
            fig_Id = fig_Id + 1;
        end
        
    end
    %     savefig(fullfile(EEGdir,sprintf('sub%03.0f_kernel_reg.fig',subID)))
    
    
end
savefig(fullfile(EEGdir,'all_subjs_kernel_reg.fig'))

%% repeat for different block conditions


for sb = 1:length(subj_list)
    subID = subj_list(sb);
    load(fullfile(EEGdir, sprintf('sub%03.0f_betas_kernel_reg.mat',subID)))
    load(fullfile(EEGdir, sprintf('sub%03.0f_block_sess_id.mat',subID)))
  
idx_1 = 1;
idx_2 = 1;
idx_3 = 1;
idx_4 = 1;

for r = 1:4
    for sess = 1:6
        for bl = 1:4
            
            ID = num2str(block_Session_ID(sess,bl,2));
            
            switch ID
                
                case '1'
                    betas_1{r}(:,:,idx_1) = squeeze(betas{r}(:,:,sess,bl));
                    idx_1 = idx_1 + 1;
                case '2'
                    
                    betas_2{r}(:,:,idx_2) = squeeze(betas{r}(:,:,sess,bl));
                    idx_2 = idx_2 + 1;
                    
                case '3'
                    
                    betas_3{r}(:,:,idx_3) = squeeze(betas{r}(:,:,sess,bl));
                    idx_3 = idx_3 + 1;
                    
                case '4'
                    
                    betas_4{r}(:,:,idx_4) = squeeze(betas{r}(:,:,sess,bl));
                    idx_4 = idx_4 + 1;
            end
        end
    end
    
betas_mean_1{sb,r} = nanmean(squeeze(betas_1{r}(channel_ind,:,:)),2);
betas_mean_2{sb,r} = nanmean(squeeze(betas_2{r}(channel_ind,:,:)),2);
betas_mean_3{sb,r} = nanmean(squeeze(betas_3{r}(channel_ind,:,:)),2);
betas_mean_4{sb,r} = nanmean(squeeze(betas_4{r}(channel_ind,:,:)),2);


end
end


for r = 1:4
    fig_id = 1;
    figure
for sb = 1:length(subj_list)
    subID = subj_list(sb);

subplot(2,4,fig_id)
fig_id = fig_id + 1;
hold on

plot(time_idx(r).timebins,betas_mean_1{sb,r} ,'g')
plot(time_idx(r).timebins,betas_mean_2{sb,r}, 'b')
plot(time_idx(r).timebins,betas_mean_3{sb,r}, 'r')
plot(time_idx(r).timebins,betas_mean_4{sb,r}, 'y')
hold off
if sb == length(subj_list)
legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL')
end
title (sprintf('%s on EEG for subj %d',time_idx(r).name, subID));
xlabel('ms')

tidyfig; 
end 
savefig(fullfile(EEGdir,sprintf('kernel_reg_%s',time_idx(r).name)))
end
%% %18
subID = 16;

load(fullfile(EEGdir, sprintf('sub%03.0f_betas_kernel_reg.mat',subID)))
fig_id = 1; 
for r = 1:4 
    
  
    
  

    for t = 1:4
        if r ~= 4
        start_time = -1000;    
        start_time = start_time + t * 500;
        else 
            start_time = -500;    
           start_time = start_time + t * 1375; 
        end 
        
        
clear bs
clear mean_b
ft_struct.time = time_idx(r).timebins;
bs = betas{r}(:,:,:,:);

% take the average across sessions
mean_b = nanmean(bs,4);
mean_b = nanmean(mean_b,3);
ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:);

cfg = [];
cfg.xlim = [start_time start_time + 500];  % time limit

% maxlim = nanmean(mean_b(:))+ 0.5 *std(mean_b(:),'omitnan');

% minlim = nanmean(mean_b(:))- 0.5 *std(mean_b(:),'omitnan');


if fig_id < 5 
    maxlim = 0.025 ; 
minlim = -0.04 ;
elseif fig_id > 4 && fig_id < 9
    maxlim = 8; 
minlim = -4 ;
elseif fig_id > 8 && fig_id < 13
    maxlim = 8; 
minlim = -4;
elseif fig_id > 12 && fig_id < 17
    maxlim = 1.6; 
minlim = -2;
end 

cfg.zlim = [minlim maxlim];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
subplot(4,4,fig_id)
fig_id = fig_id + 1;
ft_topoplotER(cfg,ft_struct); colorbar
    end

end
savefig(fullfile(EEGdir,sprintf('sub%03.0f_topo_kernel_reg',subID)))
%% for all regressors

fig_Id = 1;
figure_counter = 0;
for i = 1:length(subj_list)
    subID = subj_list(i);
    load(fullfile(EEGdir, sprintf('sub%03.0f_betas_all_reg.mat',subID)))
    
    for r = 1:6
        subplot(4,6,fig_Id)
        plotmse(squeeze(betas{r}(channel_ind,:,:)),2,time_idx(r).timebins);
        %plot(time_idx(r).timebins,squeeze(betas{r}(channel_ind,:,:)));
        
        
        xlabel(sprintf('%s on EEG',time_idx(r).name));
        tidyfig;
        if r ==1
            title(sprintf('Subject: %3.0f Channel: %s',subID, chanlbCPZ))
            
        else
            title(sprintf('Channel: %s' ,chanlbCPZ));
            
            
        end
        
        
        
        
        if fig_Id == 24
            fig_Id = 1;
            
            
            figure_counter = figure_counter + 1;
            if figure_counter == 1
                savefig(fullfile(EEGdir,'4_subjs_all_reg1.fig'))
                figure;
            else
                savefig(fullfile(EEGdir,'4_subjs_all_reg2.fig'))
                
            end
            
            
        else
            fig_Id = fig_Id + 1;
        end
        
    end
    %     savefig(fullfile(EEGdir,sprintf('sub%03.0f_all_reg.fig',subID)))
    
    
end


%%


for sb = 1:length(subj_list)
    subID = subj_list(sb);
    load(fullfile(EEGdir, sprintf('sub%03.0f_betas_all_reg.mat',subID)))
    load(fullfile(EEGdir, sprintf('sub%03.0f_block_sess_id.mat',subID)))
  
idx_1 = 1;
idx_2 = 1;
idx_3 = 1;
idx_4 = 1;

for r = 1:6
    for sess = 1:6
        for bl = 1:4
            
            ID = num2str(block_Session_ID(sess,bl,2));
            
            switch ID
                
                case '1'
                    betas_1{r}(:,:,idx_1) = squeeze(betas{r}(:,:,sess,bl));
                    idx_1 = idx_1 + 1;
                case '2'
                    
                    betas_2{r}(:,:,idx_2) = squeeze(betas{r}(:,:,sess,bl));
                    idx_2 = idx_2 + 1;
                    
                case '3'
                    
                    betas_3{r}(:,:,idx_3) = squeeze(betas{r}(:,:,sess,bl));
                    idx_3 = idx_3 + 1;
                    
                case '4'
                    
                    betas_4{r}(:,:,idx_4) = squeeze(betas{r}(:,:,sess,bl));
                    idx_4 = idx_4 + 1;
            end
        end
    end
    
betas_mean_1{sb,r} = nanmean(squeeze(betas_1{r}(channel_ind,:,:)),2);
betas_mean_2{sb,r} = nanmean(squeeze(betas_2{r}(channel_ind,:,:)),2);
betas_mean_3{sb,r} = nanmean(squeeze(betas_3{r}(channel_ind,:,:)),2);
betas_mean_4{sb,r} = nanmean(squeeze(betas_4{r}(channel_ind,:,:)),2);


end
end

for r = 1:6
    fig_id = 1;
    figure
for sb = 1:length(subj_list)
    subID = subj_list(sb);

subplot(2,4,fig_id)
fig_id = fig_id + 1;
hold on

plot(time_idx(r).timebins,betas_mean_1{sb,r} ,'g')
plot(time_idx(r).timebins,betas_mean_2{sb,r}, 'b')
plot(time_idx(r).timebins,betas_mean_3{sb,r}, 'r')
plot(time_idx(r).timebins,betas_mean_4{sb,r}, 'y')
hold off
if sb == length(subj_list)
legend('ITIS INTS', 'ITIS INTL', 'ITIL INTS', 'ITIL INTL')
end
title (sprintf('%s on EEG for subj %d',time_idx(r).name, subID));
xlabel('ms')

tidyfig; 
end 
savefig(fullfile(EEGdir,sprintf('all_reg_%s',time_idx(r).name)))
end


%% 18
subID = 19;

load(fullfile(EEGdir, sprintf('sub%03.0f_betas_all_reg.mat',subID)))
fig_id = 1; 
for r = 1:6
    
  
    
  

    for t = 1:4
        if r ~= 6
        start_time = -1000;    
        start_time = start_time + t * 500;
        else 
            start_time = -500;    
           start_time = start_time + t * 1375; 
        end 
        
        
clear bs
clear mean_b
ft_struct.time = time_idx(r).timebins;
bs = betas{r}(:,:,:,:);

% take the average across sessions
mean_b = nanmean(bs,4);
mean_b = nanmean(mean_b,3);
keyboard;
ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:);

cfg = [];
cfg.xlim = [start_time start_time + 500];  % time limit

lim = quantile(mean_b(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);


cfg.zlim = [minlim maxlim];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
subplot(6,4,fig_id)
fig_id = fig_id + 1;
ft_topoplotER(cfg,ft_struct); colorbar
    end

end
savefig(fullfile(EEGdir,sprintf('sub%03.0f_topo_all_reg',subID)))

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

cfg = [];
% cfg.xlim = [0.3 0.5];  % time limit
cfg.zlim = [-2 1];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,ft_struct); colorbar

%% 

for i = 1:length(subj_list)
    subID = subj_list(i);
clear betas
betas = load(fullfile(EEGdir, sprintf('sub%03.0f_betas_all_reg.mat',subID)));
for r = 1:6
bs = []; 
mean_b = [];
bs = betas.betas{r}(:,:,:,:);

mean_b = nanmean(bs,4);
mean_b = nanmean(mean_b,3);

bs_all{r}(:,:,i) = mean_b; 

mean_all{r} = nanmean(bs_all{r},3);
end
end
%% 
fig_id = 1; 
for r = 1:6
    
  
    
  

    for t = 1:4
        if r ~= 6
        start_time = -1000;    
        start_time = start_time + t * 500;
        else 
            start_time = -500;    
           start_time = start_time + t * 1375; 
        end 
        
        
clear bs
clear mean_b
ft_struct.time = time_idx(r).timebins;



ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_all{r}(:,:);

cfg = [];
cfg.xlim = [start_time start_time + 500];  % time limit

lim = quantile(mean_all{r}(:),[0.1 0.9]);

minlim = lim(1);
maxlim = lim(2);


cfg.zlim = [minlim maxlim];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
subplot(6,4,fig_id)
fig_id = fig_id + 1;
ft_topoplotER(cfg,ft_struct); colorbar
    end

end
savefig(fullfile(EEGdir,sprintf('sub%03.0f_topo_all_reg',subID)))

subplot(6,4,1)
title('coherence jump')
subplot(6,4,2)
title('0 0.5sec')
subplot(6,4,3)
title('0.5 1sec')
subplot(6,4,4)
title('1 1.5sec')

subplot(6,4,5)
title('coherence jump level')

subplot(6,4,9)
title('prediction error')

subplot(6,4,13)
title('button press during trial')

subplot(6,4,17)
title('button press during intertrial')

subplot(6,4,21)
title('trial start')

%% 
for i = 1:24
    
    subplot(6,4,i)
    tidyfig;
end 



