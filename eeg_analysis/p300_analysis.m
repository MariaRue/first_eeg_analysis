% analysis of P300 signal for different blocks and different coherence
% levels

% start fieldtrip  and add folders with .mat files with data structure from
% fieldtrip after pre-processing
addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data
ft_defaults % start fieldtrip

addpath('/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub001/EEG/pre-processed_data_eeg_pilot001');
%% load data

data = load('sub001_clean.mat');

data_final = data.data_final;


%% across all sessions plot P300 for differenct coherence levels across parietal electrodes

coherence = unique(data_final.trialinfo);
coherence = coherence(1:4);


for i = 1 : 4
    
    idx_coh = data_final.trialinfo == coherence(i) | data_final.trialinfo == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh;
    cfg.channel = {'P7';'P5';'P1';'PZ';'P2';'P3';'P4';'P6';'P8';'P07';'P03';'P0Z';'P04';'P08'};
    cfg.avgoverchan = 'yes';
    data_coherence{i} = ft_selectdata(cfg,data_final);
    cfg = [];
    average_ERP{i} = ft_timelockanalysis(cfg,data_coherence{i});
    
end

cfg = [];
% cfg.channel = {'P7';'P5';'P1';'PZ';'P2';'P3';'P4';'P6';'P8';'P07';'P03';'P0Z';'P04';'P08'};
cfg.baseline = [-0.4 -0.1];
cfg.baselinetype = 'absolute';
cfg.layout = 'quickcap64.mat';
cfg.showlabels = 'yes';
% figure;
% hold on
% for i = 1:4
ft_singleplotER(cfg, average_ERP{:});
% black = 75% coherence, green = 50% coherence, red = 35% coherence, blue =
% 25% coherence
% end
% hold off
legend('25','35','50','70')
%% now repeat the same but for different blocks

% get block sequence and number of trials for each block
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub001/behaviour'); %add path to folder with behavioural data and stim structures

%% load S structure from all sessions

% high freq data
S1 = load('sub001_sess017_behav.mat'); S{1} = S1.S; % disregard last trial in
% last (4th block) because right at end of session and not in EEG events - figure out how why that happened
S2 = load('sub001_sess022_behav.mat');S{2} = S2.S;
S3 = load('sub001_sess023_behav.mat');S{3} = S3.S;

%% low freq data
% high freq data
S1 = load('sub001_sess018_behav.mat'); S{1} = S1.S; % disregard last trial in
% last (4th block) because right at end of session and not in EEG events - figure out how why that happened
S2 = load('sub001_sess019_behav.mat');S{2} = S2.S;
S3 = load('sub001_sess020_behav.mat');S{3} = S3.S;
S4 = load('sub001_sess025_behav.mat');S{4} = S4.S;

%%

% genearate counters for saving trial indices in condition vectors for each
% block
countITIS_INTES = 0;
countITIS_INTEL = 0;
countITIL_INTES = 0;
countITIL_INTEL = 0;

% pre-allocate space for all 4 conditions
cond1 = zeros(3,2);
cond2 = zeros(3,2);
cond3 = zeros(3,2);
cond4 = zeros(3,2);
% coutner of start index of first trial for each block in each condition
trial = 1;


% loop through sessions

for s = 1:numel(S)
    sess = S{s}; % get session
    
    % loop through blocks of a session
    % for each block have one condition vector in which we save the trial
    % indices they correspond to in the EEG data
    blocknum = 4;
% %     for low freq 
%         if s == 4 
%             blocknum = 1;
%         else 
%             blocknum = 4;
%             
%         end
    for b = 1 : blocknum
        
        
        % get block id, 1 - 4 each number corresponding to one of the following below
        block = sess.block_ID_cells{b};
        
        % get number of trials in that block
        total_trials = max(max(sess.blocks_shuffled{b}));
        
        % switch through block ids
        switch block
            
            
            % ITIS short, INTE short
            case '1'
                %update counter for this block condition
                countITIS_INTES = countITIS_INTES+1;
                    if s == 3
                       total_trials = total_trials-1;
                    end
                cond1(countITIS_INTES,:) = [trial, trial+(total_trials-1)];
             
                % ITIS short, INTE long
            case '2'
                %update counter for this block condition
                countITIS_INTEL = countITIS_INTEL+1;
                %
                    if s == 3
                        total_trials = total_trials -1;
                    end
                
                cond2(countITIS_INTEL,:) = [trial, trial+(total_trials-1)];
                
                % ITIS long, INTE short
            case '3'
                %update counter for this block condition
                countITIL_INTES = countITIL_INTES+1;
                    if s == 1 % for high freq
                       total_trials = total_trials-1;
                    end
                
                
%                 % for low freq
%                 if s == 2
%                     
%                     total_trials = 5;
%                     
%                 end
                
                
                cond3(countITIL_INTES,:) = [trial, trial+(total_trials-1)];
                
                % ITIS long, INTE long
            case '4'
                %update counter for this block condition
                countITIL_INTEL = countITIL_INTEL+1;
                
                cond4(countITIL_INTEL,:) = [trial, trial+(total_trials-1)];
                
        end % end of switch condition
        
        % update start index for next block
        trial = trial + total_trials;
        
    end  % end loop through blocks within a session
    
end % loop through sessions

%% now sort trials in EEG data

all_cond = {cond1, cond2, cond3, cond4};
%    

for b = 1 : 4 % loop through conditions, get all trial indices and select
    % data with field trip
    idx = zeros(355,1);
    block = all_cond{b};
    
    %     for low freq 
%         if b == 3 
%             sessionnum = 4;
%         else 
            sessionnum = 3;
            
%         end
    for s = 1 : sessionnum
        idx(block(s,1):block(s,2)) = 1;
        idx = logical(idx);
    end
    
    cfg = [];
    
    cfg.trials = idx;
    cfg.channel = {'P7';'P5';'P1';'PZ';'P2';'P3';'P4';'P6';'P8';'P07';'P03';'P0Z';'P04';'P08'};
    cfg.avgoverchan = 'yes';
    data_condition{b} =  ft_selectdata(cfg,data_final);
end

cfg = [];
cfg.baseline = [-0.4 -0.1];
cfg.baselinetype = 'absolute';
cfg.layout = 'quickcap64.mat';
%
figure
ft_singleplotER(cfg,data_condition{:})
legend('ITIshortINTEshort','ITIshortINTElong','ITIlongINTEshort','ITIlongINTElong')

%%

condition_1 = rem(data_condition{1}.trialinfo, 100);

demean_cond_1 = condition_1 - mean(condition_1);

%% put data into fieldtrip 'timelock' format

cfg = [];
cfg.keeptrials = 'yes';
tl{1} = ft_timelockanalysis(cfg, data_condition{1});
B = zeros(2,2000);
for t = 1:2000
    Y = tl{1}.trial(:,1,t);
    X = demean_cond_1;
    [B(:,t), ~, st{t}] = glmfit(X,Y);
    se(:,t) = st{t}.se;
    clear st
end

