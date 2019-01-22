current_user = 'MR';

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
    scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
    EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','EEG');
    BHVdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','EEG_pilot','sub003','behaviour');
    ft_defaults
end

%% convert EEG data; downsample to 100 Hz; bandpass filter 0.1-30Hz

subID = 3;
nSess = 5; %number of sessions

cd(EEGdatadir);
for i = 1:nSess
    fname_target = fullfile(EEGdatadir,...
        sprintf('fdspmeeg_sub%03.0f_sess%03.0f_fil001.mat',subID,i));
    if exist(fname_target,'file')
        D{i} = spm_eeg_load(fname_target);
    else
        S = [];
        S.dataset = fullfile(EEGdatadir,sprintf('sub%03.0f_sess%03.0f_fil001.set',subID,i));
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

%% load in behavioural data

for i = 1:nSess
    fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
    bhv{i} = load(fname_behav);
end


%% align behavioural data with EEG data

nBlocks = 4;

for i = 1:nSess %loop over sessions
    eeg_events = D{i}.events; %events in EEG data
    eob = find([eeg_events.value]==210); %end of block trigger
    if length(eob)~=nBlocks
        error('didn''t find 4 end of blocks');
    else
        eob = [0 eob];
        for b = 1:nBlocks % loop over blocks

            %a number of triggers in S.trigger_vals haven't made it into
            %the EEG data. We now try to correct this problem, making
            %'trigger_vals_eegmatch' - a version of S.trigger_vals that
            %only contains the triggers that are also in the EEG data
            trigger_vals_behav = bhv{i}.S.trigger_vals{b};
            tvind = find(trigger_vals_behav); tvlist = trigger_vals_behav(tvind); % a list of all the trigger values in S.trigger vals
            
            eeg_events_block = eeg_events(eob(b)+1:eob(b+1)); %eeg_events just corresponding to this block
            tvlist_eeg = [eeg_events_block.value]'; % a list of all trigger values that are in the EEG data
            
            if length(tvlist)<=length(tvlist_eeg) % quick sanity check - are there more really triggers in S.tvlist? (answer = yes for pilot dataset)
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

for i = 1:nSess
    eeg_events = D{i}.events; %events in EEG data
    eob = find([eeg_events.value]==210); %end of block trigger
    eob = [0 eob];
    
    for b = 1:nBlocks
        eeg_events_block = eeg_events(eob(b)+1:eob(b+1)); %eeg_events just corresponding to this block
        blocklength = eeg_events_block(end).time-eeg_events_block(1).time; %length of block, in seconds
        nSamplesEEG = round(blocklength*D{1}.fsample)+1; %number of samples in block, at 100 Hz (fsample)
        
        %make a vector from the eeg trigger channel that corresponds to S.trigger_vals (with 0s where no trigger occurs)
        trigger_vals_eeg = zeros(nSamplesEEG,1);
        for s = 1:length(eeg_events_block)
            trigger_vals_eeg(round((eeg_events_block(s).time-eeg_events_block(1).time)*D{i}.fsample)+1) = ...
                eeg_events_block(s).value;
        end
        
        trigger_vals_behav = bhv{i}.S.trigger_vals{b}; %trigger values from behaviour
        nSamplesBehav = length(trigger_vals_behav);
        
        %now 'zero-pad' the EEG triggers, as the behavioural triggers may run over in length
        zp_size = 500; %number of samples to zero pad by
        trigger_vals_eeg_zp = [zeros(zp_size,1); trigger_vals_eeg; zeros(zp_size,1)]; %zeropadded eeg triggervalues
        nSamplesEEG_zp = length(trigger_vals_eeg_zp);
        
        for c = 1:(nSamplesEEG_zp + 1 - nSamplesBehav)
            nMatch(c) = sum(trigger_vals_behav==trigger_vals_eeg_zp(c:c+nSamplesBehav-1));
%             
%            if i == 6 && b == 2
%                keyboard
%            end
        end
        
      
       % plot(nMatch); disp(b); pause; % this reveals a clear 'spike' in every session -
        % where the triggers in the EEG data match the behavioural triggers
        %but a bit strangely, it doesn't always seem
        %to be at 500 - it is sometimes up to half a
        %second earlier - Maria to investigate?
        
        [~,best_match(i,b)] = max(nMatch);
        
        %we may be out here by one sample - I can't quite work the indexing out, but
        %it won't matter in the grand scheme of things...
        block_start_eeg_idx(i,b) = findc(D{1}.time,eeg_events_block(1).time)-zp_size+best_match(i,b);
        block_end_eeg_idx(i,b)   = block_start_eeg_idx(i,b) + nSamplesBehav - 1;
        
        %grab the relevant data
        EEGdat{i}{b} = D{i}(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
    end
end


%% build 'sliding' GLM

clear betas
nChannels = 64;
nSess = 5;
for i = 1:nSess
    
    if i == 1
        
        nBlocks = 3;
    else
        nBlocks = 4;
    end
    
    for b = 1:nBlocks
        nLags = 150; %number of lags to test (100 lags = 1s)
        
        coherence = bhv{i}.S.coherence_frame{b}; %vector of coherence levels for this block
        coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
        % then *all* dots move in same direction, i.e. coherence = 1
        
        
        mean_coherence = bhv{i}.S.mean_coherence{b}; % vector of mean coherences of this block - to figure out trial periods
        
        
        %coherence = coherence(1:1000); % for piloting, delete once complete
        coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
        coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherence 'jumps'
        
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
        
        nF = length(coherence);
%         
%         regressor_list(1).value = coherence_jump; 
%         regressor_list(1).nLagsBack = 100; 
%         regressor_list(1).nLagsForward = 150;
%         regressor_list(1).name = 'coherence_jump'; 
%         
%         regressor_list(2).value = coherence_jump_level; 
%         regressor_list(2).nLagsBack = 100; 
%         regressor_list(2).nLagsForward = 150;
%         regressor_list(2).name = 'coherence_jump_level'; 
        
        regressor_list(1).value = button_press; 
        regressor_list(1).nLagsBack = 150; 
        regressor_list(1).nLagsForward = 150;
        regressor_list(1).name = 'button_press'; 
        
        regressor_list(2).value = trial_start;
        regressor_list(2).nLagsBack = 50; 
        regressor_list(2).nLagsForward = 200;
        regressor_list(2).name = 'trial start'; 
        
        regressor_list(4).value = EEGdat{i}{b}(63,:)';
        regressor_list(4).nLagsBack = 0; 
        regressor_list(4).nLagsForward = 0;
        regressor_list(4).name = 'confound_EOG_reg_ver'; 
        
        regressor_list(3).value = EEGdat{i}{b}(64,:)';
        regressor_list(3).nLagsBack = 0; 
        regressor_list(3).nLagsForward = 0;
        regressor_list(3).name = 'confound_EOG_reg_hor'; 
        

        
        Fs = D{i}.fsample; 
        [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs); 
   
        tmp = (pinv(lagged_design_matrix')*EEGdat{i}{b}')';
        
        for r = 1:length(regressor_list)
            betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
        end
                
    end
end

channel_ind = 40; %channel of interest (CPz = 40);
chanlabel = D{1}.chanlabels(channel_ind); chanlabel = chanlabel{1};

for r = 1:3
    figure;
    plotmse(squeeze(betas{r}(channel_ind,:,:)),2,time_idx(r).timebins);
    xlabel(sprintf('Influence of %s on EEG at time (t+X) ms',time_idx(r).name));
    title(sprintf('Channel: %s',chanlabel));
    tidyfig;
end

% 
% figure;
% plotmse(squeeze(betas_coi(:,4,:)),2,time_label2);
% xlabel('Influence of button press during coherent motion on future EEG at time t');
% title(sprintf('Channel: %s',chanlabel));
% tidyfig;
% 
% figure;
% plotmse(squeeze(betas_coi(:,5,:)),2,time_label2);
% xlabel('Influence of button press during incoherent motion on future EEG at time t');
% title(sprintf('Channel: %s',chanlabel));
% tidyfig;
% 
% figure;
% plotmse(squeeze(betas_coi(:,6,:)),2,time_label2);
% xlabel('Influence of trial period on future EEG at time t');
% title(sprintf('Channel: %s',chanlabel));
% tidyfig;

%% %% open fieldtrip - make topoplot of regressors with fieldtrip 


% start fieldtrip  and add folders with .mat files with data structure from
% fieldtrip after pre-processing
bs = betas(:,:,:,:);

% take the average across sessions 
mean_b = mean(bs,4);

ft_defaults % start fieldtrip

ft_struct.dimord = 'chan_time';

ft_struct.label = D{1}.chanlabels; 
% ft_struct.elec = average_ERP{1}.elec;
ft_struct.avg = mean_b(:,:,3); 
ft_struct.time = time_label; 

%% plot topoplots 


cfg = [];                            
% cfg.xlim = [0.3 0.5];  % time limit               
% cfg.zlim = [0 6e-14];  % colour limit            
cfg.layout = 'quickcap64.mat';          
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,ft_struct); colorbar





%% repeat the above GLM analysis for different block types 

clear betas
clear betas_test
nChannels = 64;
nSess = 5;

for i = 1:nSess
    
    if i == 1
        
        nBlocks = 3;
    else
        nBlocks = 4;
    end
    
    for b = 1:nBlocks
        nLags = 150; %number of lags to test (100 lags = 1s)
        
      
        blockID(i,b) = str2num(bhv{i}.S.block_ID_cells{b});
      
        
        coherence = bhv{i}.S.coherence_frame{b}; %vector of coherence levels for this block
        coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
        % then *all* dots move in same direction, i.e. coherence = 1
        
        
        mean_coherence = bhv{i}.S.mean_coherence{b}; % vector of mean coherences of this block - to figure out trial periods
        
        
        %coherence = coherence(1:1000); % for piloting, delete once complete
        coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
        coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherence 'jumps'
        
        integration_start = abs([0; diff(mean_coherence)])>0; %vector of trial starts
        
        button_press = trigger_vals_eegmatch{i}{b} == 201 |... 
                       trigger_vals_eegmatch{i}{b} == 202;
                         % vector of button presses during trial periods
        %                        trigger_vals_eegmatch{i}{b} == 205 |...
        %                        trigger_vals_eegmatch{i}{b} == 206; 
        
        button_press_incoh_motion = trigger_vals_eegmatch{i}{b} == 205 |...
                                    trigger_vals_eegmatch{i}{b} == 206; 
                                % vector of button presses during intertrial periods
        
        
        nF = length(coherence);
        
        clear reg;
        for l = 1:nLags
            lback = l-1; %this is how many 'back' in time this regressor looks
            
            % regressor for ERP
            reg_zp = [zeros(nLags,1); coherence_jump]; %zero pad regressor
            reg(:,l) =  reg_zp(nLags-lback+1:length(reg_zp)-lback);
            
            % regressor for influence of coherence magnitude of step 
            reg_zp = [zeros(nLags,1); coherence_jump_level]; %zero pad regressor
            reg(:,l+nLags) =  reg_zp(nLags-lback+1:length(reg_zp)-lback);
            
            % regressor for influence of button press in trial period
            reg_zp = [zeros(nLags,1); button_press]; %zero pad regressor
            reg(:,l+2.*nLags) =  reg_zp(nLags-lback+1:length(reg_zp)-lback);
            
            
            % regressor for how button press going forward in trial period
            % - probably failed 
            reg_zp = [button_press; zeros(nLags,1); ]; %zero pad regressor
            reg(:,l+3.*nLags) =  reg_zp(l:length(reg_zp)-nLags + lback);
            
            
            % regressor for how button press going forward in intertrial
            % period - probably failed 
            reg_zp = [button_press_incoh_motion; zeros(nLags,1); ]; %zero pad regressor
            reg(:,l+4.*nLags) =  reg_zp(l:length(reg_zp)-nLags + lback);
            
            
            % regressor for influence of last trial period on signal -
            % failed 
            reg_zp = [mean_coherence; zeros(nLags,1); ]; %zero pad regressor
            reg(:,l+5.*nLags) =  reg_zp(l:length(reg_zp)-nLags + lback);
            
            
            nReg = 6; %number of regressors
            
            confound_EOG_reg = EEGdat{i}{b}(63:64,:)'; %EOG channels - include as confound regressor
            
        end
        
        
        
        tmp = (pinv([reg confound_EOG_reg])*EEGdat{i}{b}')';
        betas(:,:,:,i,b) = reshape(tmp(:,1:nLags*nReg),nChannels,nLags,nReg); %channels * lags * regressors * sessions * blocks
    end
end

condition{1} = 'ITIs INTs'; 
condition{2} = 'ITIs INTL'; 
condition{3} = 'ITIL INTs'; 
condition{4} = 'ITIL INTL'; 

for block = 1:4

    
 idx = blockID == block;    
betas_test = betas(:,:,:,idx);
    
    
time_label = 0:-10:-(nLags-1)*10;
time_label2 = 0:10:(nLags-1)*10;
channel_ind =40; %channel of interest (CPz = 40);
chanlabel = D{1}.chanlabels(channel_ind); chanlabel = chanlabel{1};

betas_coi = squeeze(betas_test(channel_ind,:,:,:,:));
betas_coi = betas_coi(:,:,:); %collapse across sessions/blocks

figure;
plotmse(squeeze(betas_coi(:,1,:)),2,time_label);
xlabel('Influence of coherence jump at time (t-X) ms on EEG at time t');
title(sprintf('Channel: %s Block: %s',chanlabel, condition{block}));
tidyfig;

figure;
plotmse(squeeze(betas_coi(:,2,:)),2,time_label);
xlabel('Influence of coherence jump magnitude at time (t-X) ms on EEG at time t');
title(sprintf('Channel: %s Block: %s',chanlabel, condition{block}));
tidyfig;

figure;
plotmse(squeeze(betas_coi(:,3,:)),2,time_label);
xlabel('Influence of button press at time (t-X) ms on EEG at time t');
title(sprintf('Channel: %s Block: %s',chanlabel, condition{block}));
tidyfig;

% figure;
% plotmse(squeeze(betas_coi(:,4,:)),2,time_label2);
% xlabel('Influence of button press during coherent motion on future EEG at time t');
% title(sprintf('Channel: %s Block: %s',chanlabel, condition{block}));
% tidyfig;
% 
% figure;
% plotmse(squeeze(betas_coi(:,5,:)),2,time_label2);
% xlabel('Influence of button press during incoherent motion on future EEG at time t');
% title(sprintf('Channel: %s Block: %s',chanlabel,condition{block}));
% tidyfig;
% 
% figure;
% plotmse(squeeze(betas_coi(:,6,:)),2,time_label2);
% xlabel('Influence of trial period on future EEG at time t');
% title(sprintf('Channel: %s Block: %s',chanlabel, condition{block}));
% tidyfig;

end


