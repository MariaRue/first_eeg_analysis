 function [session_data, combined_data, coherences, PMF_fit] = TTT_training_analysis(list, sep)

% Function to read in a list of *.rec datasets from T&T touchscreen
% training sessions and plot PMFs and RTs by stimulus coherence and
% experimental condition (where relevant).
%
% INPUTS
%   list : A text file list of full paths to the *.rec raw datasets (use
% UNIX "find" command if necessary).
%   sep : 0 = plot combined data only; 1 = plot individual sessions as well
% as combined data.
%
% NB. Can quickly print a list of full paths to the *.rec raw data using
% the following command: find /path/ . -name \*.rec -print > list.txt
%
% The script detects the type of task from the first *.rec input filename:
% fixed duration (FD); reaction time (RT); reward bias (REW); or social
% 2AFC task (SOCIAL).
%
% All datasets are reformatted into a "data" variable of the following
% format:
%
% Column    Information
% 1 :       Trial condition*
% 2 :       Choice (0 = leftwards choice, 1 = rightwards choice)
% 3 :       RT (seconds)
% 4 :       Coherence (%)
% [5 :       REW or SOCIAL cue condition ("right" = 1, "left" = 2)        ]
% [6 :       Social movie view angle ("right angle" = 1, "left angle" = 2)]
% [7 :       Social movie ID (integer)                                    ]
%
% Column 5 is only calculated for REW and SOCIAL tasks; Columns 6 & 7 are
% only calculated for SOCIAL tasks.
%
% * NB. There are n = n_coh*n_conds trial conditions. Trial condition
% numbers 1:n are for the "first" condition ("rightwards") and numbers
% (n+1):2n are for the "second" condition ("leftwards").
%
% * NB2. The convention is that "rightwards" choices and condition indexes
% are always ODD ("1") and "leftwards" ones are always EVEN ("0" or "2").
%
% ANALYSIS
% For FD and RT tasks, a script "TTT_plot_basic_PMFs_RTs.m" will plot the
% PMFs and RTs. For FD tasks, the RT is calculated from target onset. For
% RT tasks, the RT is calculated from stimulus onset.
%
% For REW and SOCIAL tasks, a script "TTT_plot_2conds_PMFs_RTs.m" will plot
% the PMFs and RTs separately for "left" and "right" cued conditions.
% Additionally, the script will calculate where there is any significant
% difference between PMFs in the two conditions (using chi-squared test);
% similarly for median RTs.
%
% PMF and RT plots will be generated for all the sessions in the list
% combined. If input argument "sep" = 1, then plots will be generated for
% each individual session as well (otherwise sep should be set to 0).
%
% At the end of the script there will be the option to generate plots for a
% user-specified subset of the input datasets.
%
% Nela Cicmil, 9th December 2016, University of Oxford (DPAG)

%% NB: (Maria)
% see code line 159/160 I addapted the code so that as for the MET data
% the fith column contains 0 for incorrect and 1 for correct choices. I did
% that since I didnt look at the social data. 


%%

% Dummy output values for debugging:
output = 0;


% Read in list of paths to *.rec datasets to analyse:

file_id = fopen(list, 'r'); % Get file identifiers
raw_dataset_files = textscan(file_id, '%s'); % Read file as string
n_datasets = length(raw_dataset_files{1});
fprintf('\n%d datasets to upload.\n', n_datasets); % Print number of datasets to analyse


% Discern type of dataset(1 = FD; 2 = RT; 3 = REW ; 4 = SOCIAL):
% NB. This assumes pattern '_[a-z]\.' is found only in dataset filename and not
% anywhere in the path.
test_str = raw_dataset_files{1}{1};
exp1 = '_[a-z]{2}\.';
exp2 = '[a-z]{2}';
matches = regexp( regexp(test_str, exp1, 'match') , exp2, 'match' );
data_code = char(matches{1});
switch data_code
    case 'fd'
        data_type = 1;
        fprintf('\nThis is a fixed duration (FD) task.\n');
    case 'rt'
        data_type = 2;
        fprintf('\nThis is a reaction time (RT) task.\n');
    case 'rw'
        data_type = 3;
        fprintf('\nThis is a reward bias (REW) task.\n');
    case 'so'
        data_type = 4;
        fprintf('\nThis is a social bias (SOCIAL) task.\n');
    otherwise
        fprintf('\n\Error: cannot determine data type.n');
end


% Get the session number for each dataset:

full_expression = '\.[0-9]{3}\.';
expression = '[0-9]{3}';

for i = 1:n_datasets
    session_ID_cell = regexp( regexp(raw_dataset_files{1}{i}, full_expression, 'match') , expression, 'match');
    session_ID(i) = str2double(session_ID_cell{1});
end



% Import each raw (binary) dataset:

fprintf('\nImporting the following dataset session numbers (this will take some time):\n');

for i = 1:n_datasets
    raw_data(i).r = exp_utils_read_data_binary(raw_dataset_files{1}{i}, 'alldata', 0);
    fprintf('\n%d\n', session_ID(i));
end



% Format raw data into columnar "data" variable as described above,
% separately for each training session:

for i = 1:n_datasets
    
    % Clear previously used variables:
    idx_completed = [];
    n_trials_per_block = [];
    idx_first_trial_per_block = 1;
    
    % For each block,
    n_blocks = raw_data(i).r.nblocks;
    
    for b = 1:n_blocks
        % Get the index for the completed trials:
        idx_completed(b).idx = union( find( strcmp(raw_data(i).r.resp.type(b,:), 'correct')) , find( strcmp(raw_data(i).r.resp.type(b,:), 'incorrect')) );
        
        % Get the number of completed trials in each block and index of first trial
        % of each block:
        n_trials_per_block(b) = size(idx_completed(b).idx, 2);
        if b > 1
            idx_first_trial_per_block(b) = idx_first_trial_per_block(b-1) + n_trials_per_block(b-1);
        end
    end
    idx_first_trial_per_block(n_blocks+1) = idx_first_trial_per_block(n_blocks) + n_trials_per_block(n_blocks);
    
    % For all dataset types, get choices from each block and append to the
    % data structure (Column 2); get reaction times (sec) from each block
    % and append to the data structure (Column 3); get stimulus coherence
    % values and append to the data structure (Column 4).
    for b = 1:n_blocks
        if n_trials_per_block(b) > 0
            session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 2) = ([raw_data(i).r.resp.choice{b, idx_completed(b).idx}]' +1) /2;
            session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 3) = ([raw_data(i).r.resp.reactiontime{b, idx_completed(b).idx}]' /1000); % Reaction times in seconds
            session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 4) = [raw_data(i).r.stim.stimulus.coherence{b, idx_completed(b).idx}]';
            % January 2017: 5th column displays correct and incorrect
            % choices. MR
            session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 5) = [raw_data(i).r.resp.correct{b, idx_completed(b).idx}]';
        end
    end
    
    % Fill in additional columns according to the data type:
    if data_type < 3
        % If FD or RT data, fill in Column 1 with coherence condition number:
        total_n_trials = size(session_data(i).data, 1);
        coh_range = unique(session_data(i).data(:, 4));
        session_data(i).coherences = coh_range;
        for t = 1:total_n_trials
            session_data(i).data(t, 1) = find(coh_range == session_data(i).data(t, 4));
        end
        
    elseif data_type == 3
        % If REW data, fill in Column 1 with reward*coh condition number,
        % and fill in Column 5 with rew bias condition number.
        
    elseif data_type == 4
        % If SOCIAL data, fill in Column 1 with social*coh condition
        % number, and fill in Column 5 with social bias condition number.
        % Fill in Column 6 with movie angle and Column 7 with movie cue ID
        % (see above).
        for b = 1:n_blocks
            if n_trials_per_block(b) > 0
                if session_ID(i) < 260
                    % Column 5: social bias condition: for sessions with
                    % 50/50 chance of R/L social cue across all stimulus
                    % conditions, in the raw data soc_cond = 1 means LEFT
                    % cue, and soc_cond = 2 means RIGHT cue. So we need to
                    % change that to fit our data formatting conventions,
                    % that is, "right" = 1, "left" = 2:
                    session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 5) = abs([raw_data(i).r.soc_cond{b, idx_completed(b).idx}]' -2) + 1;
                else
                    % Column 5: social bias condition: for sessions with
                    % different probabilites of R/L social cue across
                    % different stimulus conditions (user-specified), in
                    % the raw data soc_cond = 0 means LEFT cue, and
                    % soc_cond = 1 means RIGHT cue. So we need to change
                    % that to fit our data formatting conventions, that is,
                    % "right" = 1, "left" = 2:
                    session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 5) = abs([raw_data(i).r.soc_cond{b, idx_completed(b).idx}]' -1) + 1;
                end
                % Column 7: randomised movie cue ID
                session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1 , 7)  = [raw_data(i).r.rand_movie_index{b, idx_completed(b).idx}]';
                % Column 6: movie angle
                if max([raw_data(i).r.rand_movie_index{b, :}]) == 7
                    session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1, 6) = ones(n_trials_per_block(b), 1); % RHS viewing angle "1"
                else
                    session_data(i).data( idx_first_trial_per_block(b) : idx_first_trial_per_block(b+1)-1, 6) = ones(n_trials_per_block(b), 1) *2; % LHS viewing angle "2"
                end
            end
        end
        % Column 1: cue*coherence condition IDs
        % NB. Numbers 1:n are for the first, "rightwards", condition, and
        % numbers n+1:2n are for the second, "leftwards", condition.
        total_n_trials = size(session_data(i).data, 1);
        coh_range = unique(session_data(i).data(:, 4));
        session_data(i).coherences = coh_range;
        n_coh = size(coh_range, 1);
        for t = 1:total_n_trials
            session_data(i).data(t, 1) = find(coh_range == session_data(i).data(t, 4)) + n_coh*(session_data(i).data(t, 5) -1);
        end
    end
    
    
end



% Delete temporary raw data variables (to free up space):
clear raw_data;

% Concatenate datasets into a combined data variable:
% NB. For now we assume cue*coherence condition IDs are consistent across
% multiple sessions:
combined_data = vertcat(session_data(:).data);
coherences = unique(vertcat(session_data(:).coherences));


% Plot PMFs and RTs for combined data, and for individual sessions if
% selected:
data_info.combined = 1;
data_info.sessions = [session_ID(1) session_ID(end)];

if data_type <= 2
    % For FD and RT datasets:
    [PMF_fit.mu, PMF_fit.sd] = TTT_basic_plot_PMF_RTs(combined_data, coherences, data_info);
    if sep
        data_info.combined = 0;
        for s = 1:n_datasets
            data_info.session = session_ID(s);
            [PMF_fit.sep(s).mu, PMF_fit.sep(s).sd] = TTT_basic_plot_PMF_RTs(session_data(s).data, coherences, data_info);
        end
    end
    
elseif data_type >= 3
    % For REW and SOCIAL datasets:
    [logit_right, logit_left, PMF_fit.mu, PMF_fit.sd] = TTT_2conds_plot_PMFs_RTs(combined_data, coherences, data_info);
    if sep
        data_info.combined = 0;
        for s = 1:n_datasets
            data_info.session = session_ID(s);
            [session_data(s).logit_right, session_data(s).logit_left, PMF_fit.sep(s).mu, PMF_fit.sep(s).sd] = TTT_2conds_plot_PMFs_RTs(session_data(s).data, coherences, data_info);
        end
    end
end

% For REW and SOCIAL dataset types, statsitically test for significant
% shift of PMFs, using chi-squared test (assuming same SDs):
if data_type >= 3
    logit = [];
    logit.n = [logit_right.n logit_left.n];
    logit.x = [logit_right.x logit_left.x];
    logit.resps = [logit_right.resps logit_left.resps];
    logit.expno = [logit_right.expno logit_left.expno];
    
    [PMF_fit.stat_output] = run_2conds_Gauss_fitting(logit);
    title(sprintf('PMF shift analysis for sessions %d to %d. chi = %.2f; p = %.3f.', session_ID(1), session_ID(end), PMF_fit.stat_output.chi, PMF_fit.stat_output.pval), 'FontSize', 14);
    
    if sep
        for s = 1:n_datasets
            session_data(s).logit.n = [session_data(s).logit_right.n session_data(s).logit_left.n];
            session_data(s).logit.x = [session_data(s).logit_right.x session_data(s).logit_left.x];
            session_data(s).logit.resps = [session_data(s).logit_right.resps session_data(s).logit_left.resps];
            session_data(s).logit.expno = [session_data(s).logit_right.expno logit_left.expno];
            
            [session_data(s).stat_output] = run_2conds_Gauss_fitting(session_data(s).logit);
            title(sprintf('PMF shift analysis for session %d. chi = %.2f; p = %.3f.', session_ID(s), session_data(s).stat_output.chi, session_data(s).stat_output.pval), 'FontSize', 14);
        end
    end
end


fprintf('\nBye! :)\n');

end


% Plot PMFs and RTs for specific combined sessions (if selected):

%session_idx = input('\nInput array of session indices to combine and analyse \n(e.g. [a b c]), or 0 if not required --> ');

% Further statistical analysis to do? Basic DMAT model fitting?










