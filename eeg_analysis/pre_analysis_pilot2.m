% analysis for coherence jumps in second pilot 
addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data 
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/pre_process_eeg_pilot003');
ft_defaults
%% 

fid = fopen('data_files.txt'); % list of files I want to append and analyse 
f = textscan(fid,'%s'); 

% loop through files extract coherent motion events and append files to
% each other 
num_sess = size(f{:},1); 

gap = 10000;
for i = 1:num_sess
    disp(i);
    dataset = load(f{1}{i}); 
data{i} = dataset.data_trial;
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

%% mark channels in which is a nan period in a trial to exclude them 




for i = 1:length(data_trial.trial)
    
    vec_good_channels(:,i) = sum(isnan(data_trial.trial{i}(:,:)),2) == 0;
    
    
end 


vec_good_channels(62:end,:) = 0; 
%% sort data 

cor_right{1} = find(data_append.trialinfo == 24); 
cor_right{2} = find(data_append.trialinfo == 25); 
cor_right{3} = find(data_append.trialinfo == 26); 

cor_left{1} = find(data_append.trialinfo == 34); 
cor_left{2} = find(data_append.trialinfo == 35); 
cor_left{3} = find(data_append.trialinfo == 36); 


%%  ft_timelock 


% highest coherence 


% right 
for i = 1:length(cor_right{1})
    tr = cor_right{1}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40); 
    data_cor_right_high(i,:) = data_append.trial{tr}(40,:); 
   
end


% left 
for i = 1:length(cor_left{1})
    tr = cor_left{1}(i); 
    channels = zeros(65,1);
     channels(40) = channels(40) ; 
    data_cor_left_high(i,:) = data_append.trial{tr}(40,:); 
   
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
    channels(40) = channels(40) ; 
    data_cor_right_mid(i,:) = data_append.trial{tr}(40,:); 
   
end


% left 
for i = 1:length(cor_left{2})
    tr = cor_left{2}(i); 
    channels = zeros(65,1);
     
    data_cor_left_mid(i,:) = data_append.trial{tr}(40,:); 
   
end


mean_all_mid = nanmean([data_cor_right_mid;data_cor_left_mid]); 

subplot(3,1,3)
plot(mean_all_mid,'r')
ylabel('mean voltage')
xlabel('msec')
title('low coherences')

%% low coherences 


% right 
for i = 1:length(cor_right{3})
    tr = cor_right{3}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40); 
    data_cor_right_low(i,:) = data_append.trial{tr}(40,:); 
   
end


% left 
for i = 1:length(cor_left{3})
    tr = cor_left{3}(i); 
    channels = zeros(65,1);
    channels(40) = channels(40) ; 
    data_cor_left_low(i,:) = data_append.trial{tr}(40,:); 
   
end


mean_all_low = nanmean([data_cor_right_low;data_cor_left_low]); 

subplot(3,1,2)
plot(mean_all_low,'g')

ylabel('mean voltage')
xlabel('msec')
title('medium coherences')

%% look at single channels across trials 




%% 
clear data_cor_right_high
clear data_cor_left_high 
clear data_cor_left_mid
clear data_cor_right_mid
clear data_cor_right_low
clear data_cor_left_low
