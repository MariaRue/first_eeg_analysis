current_user = 'data';

tf_analysis = 0; % do time-frequency rahter than time-domain analysis?

switch current_user
    % set up spm (LH iMac)
    case 'LH'
        [hd,sd] = get_homedir; % what is this function doing?
        addpath(genpath(fullfile(hd,'matlab','hidden_from_matlab','spm12')));
        
        scriptdir = fullfile(hd,'projects','continuous_eeg_analysis','eeg_analysis');
        EEGdatadir= fullfile(sd,'projects','continuous_RDM','EEG_pilot','sub003','EEG');
        BHVdatadir= fullfile(sd,'projects','continuous_RDM','EEG_pilot','sub003','behaviour');
        
    case 'MR'
        % set up spm (MR iMac)
        addpath('/Users/maria/Documents/matlab/spm12');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
        scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');

        % EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','eyetracker_pilot','sub001','short_session_wo_EEG');
        % BHVdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','eyetracker_pilot','sub001','short_session_wo_EEG');

        EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','EEG');
        BHVdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','behaviour');
        BHVdatadir2= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','behaviour');
        STdatadir = fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','stim');
        ft_defaults
        
    case 'data'
        addpath('/Users/maria/Documents/matlab/spm12');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
        
        
        EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','sub031','EEG');
        BHVdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','sub031','behaviour');
        BHVdatadir2= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','sub031','behaviour');
        STdatadir = fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','sub031','stim');
        ft_defaults
        
    case 'eyetrig'
                addpath('/Users/maria/Documents/matlab/spm12');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
        EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub000','eyetracker_test', 'eeg_trigger');
        BHVdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub000','eyetracker_test', 'behaviour');
end

%% convert EEG data; downsample to 100 Hz; bandpass filter 0.1-30Hz

subID = 31;
nSess = 6; %number of sessions       
cd(EEGdatadir);
for i = 1:nSess
      
%     fname_target = fullfile(EEGdatadir,...
%         sprintf('fdspmeeg_sub%03.0f_sess%03.0f_fil001.mat',subID,i));

       fname_target = fullfile(EEGdatadir,...
        sprintf('fdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        
       
    if exist(fname_target,'file')
        D{i} = spm_eeg_load(fname_target);
        O{i} = spm_eeg_load(fname_target);
    else
        S = [];

        % S.dataset = fullfile(EEGdatadir,sprintf('LHtrig1.set',subID,i));

        
        
        
       % S.dataset = fullfile(EEGdatadir,sprintf('sub%03.0f_sess%03.0f_fil001.set',subID,i));
        
           S.dataset = fullfile(EEGdatadir,sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,i));
        

        S.mode = 'continuous';
        D{i} = spm_eeg_convert(S);
        
        S = [];
        S.D = D{i};
        S.fsample_new = 100;
        D{i} = spm_eeg_downsample(S);
        
        S = [];
        S.D = D{i};
        S.band = 'bandpass';
        S.freq = [0.1 30];
        D{i} = spm_eeg_filter(S);
    end
end
cd(scriptdir);

%% time-frequency decompoisiton


if tf_analysis
    for i = 1:nSess
        fname_target = fullfile(EEGdatadir,...
            sprintf('tf_fdspmeeg_sub%03.0f_sess%03.0f_fil001.mat',subID,i));
        if exist(fname_target,'file')
            D{i} = spm_eeg_load(fname_target);
        else
            S = [];
            S.D = D{i};
            S.channels = 'All';
            S.frequencies = 2:2:40;
            S.method = 'morlet';
            S.phase = 1;
            D{i} = spm_eeg_tf(S);
        end
    end
end

%% load in behavioural data

for i = 1:nSess
    fname_behav = fullfile(BHVdatadir2,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
    bhv{i} = load(fname_behav);
    fname_stim = fullfile(STdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
    stim{i} = load(fname_stim);
end


%% align behavioural data with EEG data

nBlocks = 4;
nSess = 6;
  for i = 1:nSess %loop over sessions
      
%       
% %       
%     if i == 1 
%         
%         nBlocks = 3;
%     else
%    nBlocks = 4;
%     end
    eeg_events = D{i}.events; %events in EEG data
    eob = find([eeg_events.value]==210); %end of block trigger
    sob = find([eeg_events.value]== 11);

     
    if length(eob)~=nBlocks 
        error('didn''t find 4 end of blocks');
    elseif length(sob) ~= nBlocks
      error('didn''t find 4 start of blocks');
    else
%           eob = [0 eob];
        
%         if i == 4 
%             nBlocks = 2; 
%             
%         else 
%             nBlocks = 4; 
%         end 
        for b = 1:nBlocks% loop over blocks
           
% %             
            
            %a number of triggers in S.trigger_vals haven't made it into
            %the EEG data. We now try to correct this problem, making
            %'trigger_vals_eegmatch' - a version of S.trigger_vals that
            %only contains the triggers that are also in the EEG data
            trigger_vals_behav = bhv{i}.B.trigger_vals{b};
            tvind = find(trigger_vals_behav); tvlist = trigger_vals_behav(tvind); % a list of all the trigger values in S.trigger vals
            
          % eeg_events_block = eeg_events(eob(b)+1:eob(b+1)); %eeg_events just corresponding to this block
   eeg_events_block = eeg_events(sob(b):eob(b)); %eeg_events just corresponding to this block
            tvlist_eeg = [eeg_events_block.value]'; % a list of all trigger values that are in the EEG data
            
            if length(tvlist)<=length(tvlist_eeg)
                % keyboard; % quick sanity check - are there more really triggers in S.tvlist? (answer = yes for pilot dataset)
                 error;
            elseif tvlist(end)~=210||tvlist_eeg(end)~=210 % second sanity check - is final trigger 210 in both lists? (answer is yes, i.e. 210 always found
                error;
            end
            
            %now loop backwards from end of tvlist_eeg and try to match with tvlist
            teegind = length(tvlist_eeg);
            
            tind = length(tvlist);
            keep = []; %this will be a list of all entries in tvlist_eeg that we keep
            while tind>0 & teegind>0
                if tvlist(tind)==tvlist_eeg(teegind)
                    keep(tind) = 1;
                    teegind = teegind - 1;
                else
                    keep(tind) = 0;
                end
                tind = tind-1;
                
                %these are the three places where there are oth bugs - might be worth further investigation by Maria
                if i==1&b==4&teegind==693 % at this time point a trigger has been send to the EEG recorder that does not exist or is not defined, which was number 2, correct trigger before would have been 26
                    %   keyboard;
                elseif i==6&b==2&teegind==877 % very weird trigger 19 only occurd once in the tvlist vector at 886 - so in the past from 877 in the eeg list - maybe all triggers are sort of delayed in the eeg trig list?
                    %  keyboard
                elseif i==6&b==3&teegind==345 % trigger 8 has been recorded in the eeg recording file but that trigger doesn't exist!
                    %   keyboard
                elseif i == 1 & b == 1 && teegind == 2 
                   % keyboard; 
                end
            end % while loop  % This might also explain the shifts we find in the lag between EEG and behav data? It also seems that over time the the lag between eeg triggers and behav triggers increases from 1 to 2 frames or more
            
            if length(tvlist(find(keep)))==length(tvlist_eeg) ...
                    && all(tvlist(find(keep))==tvlist_eeg) %we've matched them all!
                trigger_vals_eegmatch{i}{b} = trigger_vals_behav;
                trigger_vals_eegmatch{i}{b}(tvind(find(~keep))) = 0; %set the other ones to 0
            else
                fprintf('Haven''t matched session %0.0f, block %0.0f\n',i,b);
            end
            
            
        end
    end
end

%% now we find the eeg data corresponding to the relevant time-periods in S, check that trigger channel is well aligned, and snip out this eeg data

nBlocks = 4; 

for i = 1:nSess
    eeg_events = D{i}.events; %events in EEG data
    eob = find([eeg_events.value]==210); %end of block trigger
     sob = find([eeg_events.value]==11);  
    %eob = [0 eob];
    
         
%         if i == 1                 
%             nBlocks = 3;
%         else         
%             nBlocks = 4;
%         end
     

    for b = 1:nBlocks
        
      
         %eeg_events_block = eeg_events(eob(b)+1:eob(b+1)); %eeg_events just corresponding to this block
        eeg_events_block = eeg_events(sob(b):eob(b)); %eeg_events just corresponding to this block
        blocklength = eeg_events_block(end).time-eeg_events_block(1).time; %length of block, in seconds
        nSamplesEEG = round(blocklength*D{1}.fsample)+1; %number of samples in block, at 100 Hz (fsample)
        
      
        %make a vector from the eeg trigger channel that corresponds to S.trigger_vals (with 0s where no trigger occurs)
        trigger_vals_eeg = zeros(nSamplesEEG,1);
        for s = 1:length(eeg_events_block)
            trigger_vals_eeg(round((eeg_events_block(s).time-eeg_events_block(1).time)*D{i}.fsample)+1) = ...
                eeg_events_block(s).value;
         
        end
                    
        
        
        trigger_vals_behav = bhv{i}.B.trigger_vals{b}; %trigger values from behaviour
        nSamplesBehav = length(trigger_vals_behav);
        
        %now 'zero-pad' the EEG triggers, as the behavioural triggers may run over in length
        zp_size = 500; %number of samples to zero pad by
        trigger_vals_eeg_zp = [zeros(zp_size,1); trigger_vals_eeg; zeros(zp_size,1)]; %zeropadded eeg triggervalues
        nSamplesEEG_zp = length(trigger_vals_eeg_zp);
        
%         if i == 2 && b == 3
%             keyboard; 
%         end 
            
       
        
        for c = 1:(nSamplesEEG_zp + 1 - nSamplesBehav)
            nMatch(c) = sum(trigger_vals_behav==trigger_vals_eeg_zp(c:c+nSamplesBehav-1));
            
            
            %
            %            if i == 6 && b == 2
            %                keyboard
            %            end
        end    
        
   
        
%         plot(nMatch); disp(b); pause; % this reveals a clear 'spike' in every session -
        % where the triggers in the EEG data match the behavioural triggers
        %but a bit strangely, it doesn't always seem
        %to be at 500 - it is sometimes up to half a
        %second earlier - Maria to investigate?
        
        [~,best_match(i,b)] = max(nMatch);
        
        %we may be out here by one sample - I can't quite work the indexing out, but
        %it won't matter in the grand scheme of things...
         block_start_eeg_idx(i,b) = findc(D{1}.time,eeg_events_block(1).time)-zp_size+best_match(i,b);
%                
        block_end_eeg_idx(i,b)   = block_start_eeg_idx(i,b) + nSamplesBehav - 1;

% block_start_eeg_idx(i,b) = findc(D{1}.time,eeg_events_block(1).time);
% block_end_eeg_idx(i,b)   = findc(D{1}.time,eeg_events_block(end).time);
           
        %grab the relevant data
%         if tf_analysis
%         EEGdat{i}{b} = D{i}(:,:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1); 
%         EEG_time{i}{b} = O{i}(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
%         else 
         EEGdat{i}{b} = D{i}(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
     %   EEGdat{i}{b} = D{i}(:,eeg_events(sob(b)).time:eeg_events(eob(b)).time,1);
%         end 
    end
end
%% compute the LRP in time domain 
% response to actual trial periods 
i = 2; 
b = 1; 
LH_trig = 205;
RH_trig = 201;

% find sample id for these responses 
LH_resp = find(trigger_vals_eegmatch{i}{b} == LH_trig); 
RH_resp = find(trigger_vals_eegmatch{i}{b} == RH_trig);



% figure out what coherence level that responds to 
for l = 1:length(LH_resp)
idx = l; 


    while bhv{i}.S.mean_coherence_org{b}(LH_resp(idx)) == 0 
        
        idx = idx - 1; 
        
    end 
    LH_coh(l) = bhv{i}.S.mean_coherence_org{b}(LH_resp(idx));
    
    eegC3L(l,:) = EEGdat{i}{b}(29,LH_resp(l) - 10 : LH_resp(l) + 150);
    eegC4L(l,:) = EEGdat{i}{b}(33,LH_resp(l) - 10 : LH_resp(l) + 150);
   
end



% figure out what coherence level that responds to 
for l = 1:length(RH_resp)
idx = l; 


    while bhv{i}.S.mean_coherence_org{b}(RH_resp(idx)) == 0 
        
        idx = idx - 1; 
        
    end 
    RH_coh(l) = bhv{i}.S.mean_coherence_org{b}(RH_resp(idx));
    
    eegC3R(l,:) = EEGdat{i}{b}(29,RH_resp(l) - 10 : RH_resp(l) + 150);
    eegC4R(l,:) = EEGdat{i}{b}(33,RH_resp(l) - 10 : RH_resp(l) + 150);
    
end
% according to Luck book C3 (left) and C4 (right) electrodes might be good
% candidates to calculate the LRP (chan 29 + chan 33)

%%

coh_list = [0.3 0.4 0.5]; 
for c = 1:3 
    
    clear ErightRleft EleftRleft EleftRright ErightRright
    clear LH_coh_trig RH_coh_trig
    
    LH_coh_trig = abs(LH_coh) == coh_list(c); 
    RH_coh_trig = abs(RH_coh) == coh_list(c); 
    
    
    ErightRleft = mean(eegC4L(LH_coh_trig,:)); 
    EleftRleft = mean(eegC3L(LH_coh_trig,:)); 
    
    EleftRright = mean(eegC3R(RH_coh_trig,:)); 
    ErightRright = mean(eegC4R(RH_coh_trig,:)); 
    
   
    LRP(c,:) = abs((ErightRleft - EleftRleft)+(EleftRright - ErightRright))/2;
    

    
end 







%% build 'sliding' GLM

clear betas
nChannels = 64;
nSess = 6;
for i = 1:nSess
           
    disp(i); 

        if i == 7
        
        nBlocks = 3;  
             blocks = [1,2,4];
    elseif i == 8
         nBlocks = 3;
         blocks = [2,3,4];
    else
        nBlocks = 4;
        blocks = [1 2 3 4];
    end
%     
    for blockcount = 1:nBlocks
        b = blocks(blockcount);
        
%         if i == 1     
%             b = b+ 1;
%         end 
         blockID(i,b) = str2num(stim{i}.S.block_ID_cells{b});
        disp(b); 
        nLags = 150; %number of lags to test (100 lags = 1s)
        
        coherence = bhv{i}.B.coherence_frame{b}; %vector of coherence levels for this block
        coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
        % then *all* dots move in same direction, i.e. coherence = 1
        
        coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
        coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherence 'jumps'
        
        mean_coherence = bhv{i}.B.mean_coherence{b}; % vector of mean coherences of this block - to figure out trial periods
        
%         
% % difference between coherence at t and t-1 (0 at the start because
% % coherence is undefined at t0 so cannot calculate diff between t0 and t1)
% coherence_differences = [0; diff(coherence)];
% 
% % absolute value of this tell us the magnitude of the jump at this time point
% coherence_jump_level = abs(coherence_differences);
% 
% % all absolute changes are positive, so >0 gives us ?did a jump occur??
% coherence_jump = coherence_jump_level > 0;
        
        integration_start = abs([0; diff(mean_coherence)])>0; %vector of trial starts
        
        button_press = trigger_vals_eegmatch{i}{b} == 201 |... % vector of button presses during trial periods
            trigger_vals_eegmatch{i}{b} == 202;
        %                        trigger_vals_eegmatch{i}{b} == 205 |...
        %                        trigger_vals_eegmatch{i}{b} == 206;
        
        button_press_incoh_motion = trigger_vals_eegmatch{i}{b} == 205 |...
            trigger_vals_eegmatch{i}{b} == 206; % vector of button presses during intertrial periods
        
        
        trial_start = trigger_vals_eegmatch{i}{b} == 30 |...  % get start of each trial for all coherence levels
            trigger_vals_eegmatch{i}{b} == 40 |...
            trigger_vals_eegmatch{i}{b} == 50 |...
            trigger_vals_eegmatch{i}{b} == 130 |...
            trigger_vals_eegmatch{i}{b} == 140 |...
            trigger_vals_eegmatch{i}{b} == 150;
        
        
        % regressor for prediciton error
        coherences = [];
        coherences = bhv{i}.B.coherence_frame{b};
        diff_coherences = diff(coherences(coherence_jump));
        diff_coherences = [coherences(1); diff_coherences]; % differnce to prev cohernce for first coherence is that coherence itself
        jump_idx = find(coherence_jump);
        coherence_level_difference = zeros(size(coherences,1),1);
        coherence_level_difference(jump_idx) = abs(diff_coherences);
        
        nF = length(coherence);
        %
        regressor_list(1).value = coherence_jump;
        regressor_list(1).nLagsBack = 100;
        regressor_list(1).nLagsForward = 150;
        regressor_list(1).name = 'coherence_jump';
        
        regressor_list(2).value = coherence_jump_level;
        regressor_list(2).nLagsBack = 100;
        regressor_list(2).nLagsForward = 150;
        regressor_list(2).name = 'coherence_jump_level';
        
        regressor_list(3).value = coherence_level_difference;
        regressor_list(3).nLagsBack = 150;
        regressor_list(3).nLagsForward = 150;
        regressor_list(3).name = 'prediction error';
        
        
        regressor_list(4).value = abs(coherences);
        regressor_list(4).nLagsBack = 100;
        regressor_list(4).nLagsForward = 150;
        regressor_list(4).name = 'abs stimulus';
        
        regressor_list(5).value = button_press;
        regressor_list(5).nLagsBack = 150;
        regressor_list(5).nLagsForward = 150;
        regressor_list(5).name = 'button_press';
        
        regressor_list(6).value = button_press_incoh_motion;
        regressor_list(6).nLagsBack = 150;
        regressor_list(6).nLagsForward = 150;
        regressor_list(6).name = 'iti button press';
        
        regressor_list(7).value = trial_start;
        regressor_list(7).nLagsBack = 50;
        regressor_list(7).nLagsForward = 500;
        regressor_list(7).name = 'trial start';
        
        regressor_list(8).value = EEGdat{i}{b}(63,:,:)';
        regressor_list(8).nLagsBack = 0;
        regressor_list(8).nLagsForward = 0;
        regressor_list(8).name = 'confound_EOG_reg_ver';
        
        regressor_list(9).value = EEGdat{i}{b}(64,:,:)';
        regressor_list(9).nLagsBack = 0;
        regressor_list(9).nLagsForward = 0;
        regressor_list(9).name = 'confound_EOG_reg_hor';
%         
        regressor_list(1:3) = [];

     
        Fs = D{i}.fsample;
        [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs);
      
        if tf_analysis 
         for freq = 1:size(EEGdat{i}{b},2)   
                    tmp = (geninv(lagged_design_matrix')*squeeze(EEGdat{i}{b}(:,freq,:))')';
        
        for r = 1:length(regressor_list)
            betas{freq}{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
        end
        end 
            
        else
        tmp = (geninv(lagged_design_matrix')*EEGdat{i}{b}')';
        
        for r = 1:length(regressor_list)
            betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
        end
        end
    end
end

savefile = fullfile(EEGdatadir,'betas_tf.mat'); 
save(savefile, 'betas'); 

channel_ind = 40; %channel of interest (CPz = 40);
chanlabel = D{1}.chanlabels(channel_ind); chanlabel = chanlabel{1};


save_name = sprintf('betas_sub%03.0f.mat',subID);
save(fullfile(EEGdatadir,save_name),'betas', 'time_idx','chanlabel','channel_ind');

if tf_analysis
frequency = 10; 
freqlabel = D{1}.frequencies(frequency);
end 
%% 


for r = 1:6
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

% plot frequency responses for each regressors as map 
if tf_analysis 
for r = 1:6 
    
    
    
    
    
    
end 
end 

%% %% open fieldtrip - make topoplot of regressors with fieldtrip


% start fieldtrip  and add folders with .mat files with data structure from
% fieldtrip after pre-processing

clear bs 
clear mean_b
ft_struct.time = time_idx(6).timebins;
bs = betas{6}(:,:,:,:);

% take the average across sessions
mean_b = mean(bs,4);
mean_b = mean(mean_b,3);
ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = D{1}.chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:);

%% plot topoplots


cfg = [];
% cfg.xlim = [0.3 0.5];  % time limit
cfg.zlim = [-2 1];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,ft_struct); colorbar





%% repeat the above GLM analysis for different block types


clear betas_test

condition{1} = 'ITIs INTs';
condition{2} = 'ITIs INTL';
condition{3} = 'ITIL INTs';
condition{4} = 'ITIL INTL';

for i = 1:4
    for n = 1:8
        
        if i == 2
        betas_test{i}{n} = nan(64,length(time_idx(n).dm_row_idx),3);
        else
        betas_test{i}{n} = nan(64,length(time_idx(n).dm_row_idx),4);
        end 
    end
    
    cond{i} = 0; 
end

for i = 1:nSess
    
    if i == 1 
        nBlocks = 3;
    else
        nBlocks = 4; 
    end 
for block = 1:nBlocks
    idx = blockID(i,block);
    

    cond{idx} = cond{idx} + 1; 
    
    
    
    for r = 1:8
    
    betas_test{idx}{r}(:,:,cond{idx}) = betas{r}(:,:,i,block);
    end 
end 
end 
    %% 
    for c = 1:4
for r = 1:6
    figure;
  
    plotmse(squeeze(betas_test{3}{r}(channel_ind,:,:)),2,time_idx(r).timebins);
  
    title(sprintf('Channel: %s' ,chanlabel));
   
    
  xlabel(sprintf('Influence of %s on EEG at time (t+X) ms',time_idx(r).name));
    tidyfig;
    
    data = betas_test{c}{r}(channel_ind,:,:); 
    cond{c}{r} = mean(squeeze(data),2); 
end
    end
    
    %%  hold on 
    hold on 
    for c = 1:4
    plot(cond{c}{6})
    legend('ITIS_INTES', 'ITIS_INTEL', 'ITIL_INTES', 'ITIL_INTEL')
    end
    
%% blot betas in topoplot for different conditions 


% start fieldtrip  and add folders with .mat files with data structure from
% fieldtrip after pre-processing
bs = betas_test{1}{3}(:,:,:,:);

% take the average across sessions
mean_b = nanmean(bs,3);

ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = D{1}.chanlabels;
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:);
ft_struct.time = time_idx(3).timebins;
%% plot topoplots for different conditions


cfg = [];
% cfg.xlim = [0.3 0.5];  % time limit
% cfg.zlim = [0 6e-14];  % colour limit
cfg.layout = 'quickcap64.mat';
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,ft_struct); colorbar





