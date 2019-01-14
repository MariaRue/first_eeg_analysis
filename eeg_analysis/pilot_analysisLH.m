%% set up spm

[hd,sd] = get_homedir;
addpath(genpath(fullfile(hd,'matlab','hidden_from_matlab','spm12')));

scriptdir = fullfile(hd,'projects','continuous_RDM','pilot_analysis');
EEGdatadir= fullfile(sd,'projects','continuous_RDM','EEG_pilot','sub003','EEG');
BHVdatadir= fullfile(sd,'projects','continuous_RDM','EEG_pilot','sub003','behaviour');


%% convert EEG data; downsample to 100 Hz; bandpass filter 0.1-30Hz

subID = 3;
nSess = 7; %number of sessions

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
                if i==1&b==4&teegind==693
                    %keyboard;
                elseif i==6&b==2&teegind==877
                    %keyboard
                elseif i==6&b==3&teegind==345
                   %keyboard 
                end
            end
            
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
        end
        %plot(nMatch);pause; % this reveals a clear 'spike' in every session -
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
for i = 1:nSess
    for b = 1:nBlocks
        nLags = 150; %number of lags to test (100 lags = 1s)
        
        coherence = bhv{i}.S.coherence_frame{b}; %vector of coherence levels for this block
        coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1 
                                                                  % then *all* dots move in same direction, i.e. coherence = 1
        
        %coherence = coherence(1:1000); % for piloting, delete once complete
        coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
        coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherence 'jumps'
        
        nF = length(coherence);
        
        clear reg;
        for l = 1:nLags
            lback = l-1; %this is how many 'back' in time this regressor looks
            
            reg_zp = [zeros(nLags,1); coherence_jump]; %zero pad regressor
            reg(:,l) =  reg_zp(nLags-lback+1:length(reg_zp)-lback);
            
            reg_zp = [zeros(nLags,1); coherence_jump_level]; %zero pad regressor
            reg(:,l+nLags) =  reg_zp(nLags-lback+1:length(reg_zp)-lback);
            
            nReg = 2; %number of regressors
            
            confound_EOG_reg = EEGdat{i}{b}(63:64,:)'; %EOG channels - include as confound regressor
        end
        
        tmp = (pinv([reg confound_EOG_reg])*EEGdat{i}{b}')';
        betas(:,:,:,i,b) = reshape(tmp(:,1:nLags*nReg),nChannels,nLags,nReg); %channels * lags * regressors * sessions * blocks
    end
end

time_label = 0:-10:-(nLags-1)*10;
channel_ind = 40; %channel of interest (CPz = 40);
chanlabel = D{1}.chanlabels(channel_ind); chanlabel = chanlabel{1};

betas_coi = squeeze(betas(channel_ind,:,:,:,:));
betas_coi = betas_coi(:,:,:); %collapse across sessions/blocks

figure;
plotmse(squeeze(betas_coi(:,1,:)),2,time_label);
xlabel('Influence of coherence jump at time (t-X) ms on EEG at time t');
title(sprintf('Channel: %s',chanlabel));
tidyfig;

figure;
plotmse(squeeze(betas_coi(:,2,:)),2,time_label);
xlabel('Influence of coherence jump magnitude at time (t-X) ms on EEG at time t');
title(sprintf('Channel: %s',chanlabel));
tidyfig;

