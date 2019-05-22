% This is a analysis pipeline for EEG data that is in SPM format but not
% preprocessed yet. This scirpt should automate the process of reading in
% the spm data and running several GLMs on the data, plotting the betas and
% topoplots

%% specify the datasets to read in by identifying the subid
% in terminal go to the data folder containing sufolders for each
% participant --> ls > txt (make sure that list contains only direct
% subject directories)
scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
addpath('/Users/maria/Documents/matlab/spm12');
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
subj_list = [16,18:21,24,26,32];

%% loop through each subject and load the data and match the behavioural triggers with the EEG triggers and save the resulting eegdatasets
for sj = 1:length(subj_list)

    clearvars -EXCEPT scriptdir EEGdir subj_list sj
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    
    
    
    cd (EEGdatadir)
    
    eeg_files = dir('*set');
    
    % if raw data hasn't been transformed into eeglab set files then run
    % the following code and save raw data as set files
    if isempty({eeg_files.name})
        addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
        eeglab
        
        eeg_raw = dir('*dpa');
        nSess = length({eeg_raw.name});
        for l = 1:nSess
            
            
            fname_load = fullfile(EEGdatadir,...
                sprintf('sub%03.0f_sess%03.0f_eeg.cdt',subID,l));
            EEG = loadcurry(fname_load, 'CurryLocations', 'False');
            
            pop_saveset(EEG,'filename',sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,l));
        end
        
    else
        
        nSess = length({eeg_files.name});
        
    end
    
    
    
    
    % load in data from one subject or created SPM file if that hasn't
    % happened yet.
    for i = 1:nSess
        
        
        
        fname_target = fullfile(EEGdatadir,...
            sprintf('fdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        
        
        if exist(fname_target,'file')
            D{i} = spm_eeg_load(fname_target);
            O{i} = spm_eeg_load(fname_target);
        else
            S = [];
            
            
            
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
        
        % load in behavioural data
        fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
        bhv{i} = load(fname_behav);
        fname_stim = fullfile(STdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
        stim{i} = load(fname_stim);
        
    end
    
    
    cd (scriptdir)
    
    
    
    % match eeg triggers and behavioural triggers if EEGdat doesn't exist.
    % EEGdat is EEG data that matches the trigger vals in the behaviour
    eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
    
    if exist(eegdat_fname) ~= 2
        
        
        nBlocks = 4;
        
        % initiliaise matrix with logicals for which blocks could be matched
        % for a session and blockID (1st - 4th Column =
        % block logical, sheet 2 = 1-4th column
        % BlockID corresponding to block in first sheet)
        % rows indicate session
        
        block_Session_ID = zeros(6,4,2);
        
        
        for i = 1:nSess %loop over sessions
            
            
            eeg_events = D{i}.events; %events in EEG data
            eob = find([eeg_events.value]==210); %end of block trigger
            sob = find([eeg_events.value]== 11);
            
            
            if length(eob)~=nBlocks
                error('didn''t find 4 end of blocks');
            elseif length(sob) ~= nBlocks
                error('didn''t find 4 start of blocks');
            else
                
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
                        
                        %                     %these are the three places where there are oth bugs - might be worth further investigation by Maria
                        %                     if i==1&b==4&teegind==693 % at this time point a trigger has been send to the EEG recorder that does not exist or is not defined, which was number 2, correct trigger before would have been 26
                        %                         %   keyboard;
                        %                     elseif i==6&b==2&teegind==877 % very weird trigger 19 only occurd once in the tvlist vector at 886 - so in the past from 877 in the eeg list - maybe all triggers are sort of delayed in the eeg trig list?
                        %                         %  keyboard
                        %                     elseif i==6&b==3&teegind==345 % trigger 8 has been recorded in the eeg recording file but that trigger doesn't exist!
                        %                         %   keyboard
                        %                     elseif i == 1 & b == 1 && teegind == 2
                        %                         % keyboard;
                        %                     end
                    end % while loop  % This might also explain the shifts we find in the lag between EEG and behav data? It also seems that over time the the lag between eeg triggers and behav triggers increases from 1 to 2 frames or more
                    
                    if length(tvlist(find(keep)))==length(tvlist_eeg) ...
                            && all(tvlist(find(keep))==tvlist_eeg) %we've matched them all!
                        trigger_vals_eegmatch{i}{b} = trigger_vals_behav;
                        trigger_vals_eegmatch{i}{b}(tvind(find(~keep))) = 0; %set the other ones to 0
                        
                        block_Session_ID(i,b,1) = 1;
                    else
                        fprintf('Haven''t matched session %0.0f, block %0.0f\n',i,b);
                        
                    end
                    
                    block_Session_ID(i,b,2) = str2double(stim{i}.S.block_ID_cells{b});
                    
                end
            end
        end
        
        save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_block_sess_id']),'block_Session_ID','trigger_vals_eegmatch');
        save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']),'trigger_vals_eegmatch');
        
        % now we find the eeg data corresponding to the relevant time-periods in S, check that trigger channel is well aligned, and snip out this eeg data
        
        
        
        for i = 1:nSess
            eeg_events = D{i}.events; %events in EEG data
            eob = find([eeg_events.value]==210); %end of block trigger
            sob = find([eeg_events.value]==11);
            
            
            for b = 1:nBlocks
                
                if block_Session_ID(i,b,1)
                    
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
                    
                    
                    
                    for c = 1:(nSamplesEEG_zp + 1 - nSamplesBehav)
                        nMatch(c) = sum(trigger_vals_behav==trigger_vals_eeg_zp(c:c+nSamplesBehav-1));
                        
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
                    
                    %
                    EEGdat{i}{b} = D{i}(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
                end
            end
        end
        save(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_EEGdat']),'EEGdat');
        
    else
        
        load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_EEGdat']));
        load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_block_sess_id']));
        load(fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']));
    end
    
    
    % run GLM


    nBlocks = 4;
    nChannels = 64;
    save_name = sprintf('sub%03.0f_betas_all_reg.mat',subID);

    if exist(fullfile(EEGdir,'preprocessed_EEG_dat',save_name)) ~= 2
        
    
    for i = 1:nSess
        
        disp(i);
        
        %
        for b = 1:nBlocks
            
            if block_Session_ID(i,b,1)
                
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
                
                regressor_list(4).value = button_press;
                regressor_list(4).nLagsBack = 150;
                regressor_list(4).nLagsForward = 150;
                regressor_list(4).name = 'button_press';
                
                regressor_list(5).value = button_press_incoh_motion;
                regressor_list(5).nLagsBack = 150;
                regressor_list(5).nLagsForward = 150;
                regressor_list(5).name = 'iti button press';
                
                regressor_list(6).value = trial_start;
                regressor_list(6).nLagsBack = 50;
                regressor_list(6).nLagsForward = 500;
                regressor_list(6).name = 'trial start';
                
                regressor_list(7).value = EEGdat{i}{b}(63,:,:)';
                regressor_list(7).nLagsBack = 0;
                regressor_list(7).nLagsForward = 0;
                regressor_list(7).name = 'confound_EOG_reg_ver';
                
                regressor_list(8).value = EEGdat{i}{b}(64,:,:)';
                regressor_list(8).nLagsBack = 0;
                regressor_list(8).nLagsForward = 0;
                regressor_list(8).name = 'confound_EOG_reg_hor';
                %
                
                
                
                Fs = D{i}.fsample;
                [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs);
                
                
                tmp = (geninv(lagged_design_matrix')*EEGdat{i}{b}')';
                
                for r = 1:length(regressor_list)
                    betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
                end
                
            end
        end
    end
    
    
    channel_ind = 40; %channel of interest (CPz = 40);
    chanlbCPZ = D{1}.chanlabels(channel_ind); chanlbCPZ = chanlbCPZ{1};
    chanlabels = D{1}.chanlabels;
    
    
  
    save(fullfile(EEGdir,'preprocessed_EEG_dat',save_name),'betas', 'time_idx','chanlbCPZ','channel_ind','chanlabels');
    end
    clear betas regressor_list

    for i = 1:nSess
        
        disp(i);
        
        %
        for b = 1:nBlocks
            disp(b)
            if block_Session_ID(i,b,1)
                
                
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

                
    
                %
                
                regressor_list(1).value = abs(coherences);
                regressor_list(1).nLagsBack = 100;
                regressor_list(1).nLagsForward = 150;
                regressor_list(1).name = 'abs stimulus';
                
                regressor_list(2).value = button_press;
                regressor_list(2).nLagsBack = 150;
                regressor_list(2).nLagsForward = 150;
                regressor_list(2).name = 'button_press';
                
                regressor_list(3).value = button_press_incoh_motion;
                regressor_list(3).nLagsBack = 150;
                regressor_list(3).nLagsForward = 150;
                regressor_list(3).name = 'iti button press';
                
                regressor_list(4).value = trial_start;
                regressor_list(4).nLagsBack = 50;
                regressor_list(4).nLagsForward = 500;
                regressor_list(4).name = 'trial start';
                
                regressor_list(5).value = EEGdat{i}{b}(63,:,:)';
                regressor_list(5).nLagsBack = 0;
                regressor_list(5).nLagsForward = 0;
                regressor_list(5).name = 'confound_EOG_reg_ver';
                
                regressor_list(6).value = EEGdat{i}{b}(64,:,:)';
                regressor_list(6).nLagsBack = 0;
                regressor_list(6).nLagsForward = 0;
                regressor_list(6).name = 'confound_EOG_reg_hor';
                
                
                
                 Fs = D{i}.fsample;
                [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs);
                
                
                tmp = (geninv(lagged_design_matrix')*EEGdat{i}{b}')';
                
                for r = 1:length(regressor_list)
                    betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
                end
                
            end
        end
    end
    
    channel_ind = 40; %channel of interest (CPz = 40);
    chanlbCPZ = D{1}.chanlabels(channel_ind); chanlbCPZ = chanlbCPZ{1};
    chanlabels = D{1}.chanlabels;
    
    
    save_name = sprintf('sub%03.0f_betas_kernel_reg.mat',subID);
    save(fullfile(EEGdir,'preprocessed_EEG_dat',save_name),'betas', 'time_idx','chanlbCPZ','channel_ind','chanlabels');

    
end % loop through subject folders


