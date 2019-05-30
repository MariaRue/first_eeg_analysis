function [raw_data, formatted_data] = social_pilot_analysis(list)

% Function to read in a list of *.rec datasets and plot PMFs and RTs for
% two social conditions, for individual datasets and all datasets combined.
%
% Matlab datasets have the following format:
%
% condition_ID   resp   RT(ms)   in/correct   coherence   social_cond   movie_angle   cue_ID
%
% where:
% - condition_ID : number given according to social condition and coherence
% value
% - resp(0 or 1) : 1 for Rightwards choice, 0 for Leftwards choice
% - RT : reaction time in seconds
% - in/correct: 1 for "correct" (rewarded), 0 for "incorrect" (unrewarded)
% - coherence : motion coherence of stimulus dots, -ve = Leftward motion,
% +ve = Rightwards motion
% - social_cond : 1 for Rightwards social cue, 0 for Leftwards social cue
% - movie_angle : 11 for RHS shoulder angle, 10 for LHS shoulder angle
% - cue_ID : the ID number for the list of movies for the given movie angle
%
% Also option to plot data separately for each video cue (for data combined
% over all fiels using the same video cue set).
%
%
% Nela Cicmil, 14th November 2016, University of Oxford


% Read in the list of datasets to analyse:

file_id = fopen(list, 'r'); % Get file identifier
dataset_files = textscan(file_id, '%s'); % Read file as string
n_datasets = length(dataset_files{1});
fprintf('\n%d datasets to upload:', n_datasets);


% Import each raw (binary) dataset:

for n = 1:n_datasets
    
    fprintf('\n%d', n);
    fprintf('\t%s', dataset_files{1}{n});
    raw_data(n).r = exp_utils_read_data_binary(dataset_files{1}{n}, 'alldata', 0);
    
end

fprintf('\n');


% Format each raw dataset into a data structure (as described in function info text above):
fprintf('\n%d datasets to analyse:', n_datasets);
for n = 1:n_datasets
    fprintf('\n%d', n);
    fprintf('\t%s', dataset_files{1}{n});
    
    r = [];
    r = raw_data(n).r;
    
    % Clear preivously used variables:
    idx_completed = [];
    n_trials_pb = [];
    idx_first_trial_pb = [];
    data = [];
    
    % For each block, get index for the *completed* trials:
    for b = 1:r.nblocks
        idx_completed(b).idx = union( find(strcmp(r.resp.type(b,:),'correct')) , find(strcmp(r.resp.type(b,:),'incorrect')) );
        %idx_completed(b).idx = find(abs([r.resp.choice{b,:}]) == 1);
    end
    
    % Get the number of trials in each block and index of first trial of
    % each block:
    idx_first_trial_pb = 1;
    for b = 1:r.nblocks
        n_trials_pb(b) = size(idx_completed(b).idx,2);
        if b > 1
            idx_first_trial_pb(b) = idx_first_trial_pb(b-1) + n_trials_pb(b-1);
        end
    end
    idx_first_trial_pb(r.nblocks+1) = idx_first_trial_pb(r.nblocks) + n_trials_pb(r.nblocks);
    
    % Get choices from each block and append to the data structure (Column 2):    
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            %all_resps_temp = [r.resp.choice{b,:}];
            %resps_temp = all_resps_temp(idx_completed(b).idx);
            %data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 2) = resps_temp';
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 2) = ([r.resp.choice{b, idx_completed(b).idx}]' +1) / 2 ;
        end
    end
    
    % Get reaction times from each block and append to the data structure (Column 3):
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            %RTs_temp = [r.resp.reactiontime{b,:}]; % Reaction times in milliseconds (ms)
            %data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 3) = RTs_temp';
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 3) = [r.resp.reactiontime{b, idx_completed(b).idx}]' ; % still in ms (raw format units)
        end
    end
    
    
    % Get whether trial was correct (rewarded) or incorrect (unrewarded)
    % (Column 4):
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 4) = [strcmp(r.resp.type(b, idx_completed(b).idx), 'correct')]' ;
        end
    end
    
    % Get coherence values from each block and append to the data
    % structure (Column 5):
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            %all_cohvals_temp = [r.stim.stimulus.coherence{b,:}];
            %cohvals_temp = all_cohvals_temp(idx_completed(b).idx);
            %data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 4) = cohvals_temp';
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 5) = [r.stim.stimulus.coherence{b, idx_completed(b).idx}]' ;
        end
    end
    
    
    % Get trial social conditions from each block and append to the data
    % structure (Column 6):
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            %all_soc_conds_temp = [r.soc_cond{b,:}];
            %soc_conds_temp = all_soc_conds_temp(idx_completed(b).idx);
            %data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 5) = soc_conds_temp';
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 6) = [r.soc_cond{b, idx_completed(b).idx}]'-1 ;
        end
    end
    
    
    % Get movie angle for the dataset and append to the data structure (Column 7):
    if max([r.rand_movie_index{1,:}]) > 5 % If RHS angle
        data(:,7) = ones(size(data,1), 1) * 11;
    else
        data(:,7) = ones(size(data,1), 1) * 10; % if LHS angle (fewer movies)
    end
    
    
    % Get trial movie_ID from each block and append to the data structure (Column 8):
    for b = 1:r.nblocks
        if n_trials_pb(b) > 0
            data( idx_first_trial_pb(b) : idx_first_trial_pb(b+1)-1 , 8) = [r.rand_movie_index{b, idx_completed(b).idx}]' ;
        end
    end
    
    formatted_data(n).data = data;
    

    
end
fprintf('\n');

% If multiple datasets, combine into a combined dataset:
if n_datasets > 2
    
end



% Plotting PMF for each dataset: 
for n = 1:n_datasets
    
    data = formatted_data(n).data;
    
    % Separate data according to social condition: 
    idx_leftcue = find(data(:,6) == 0); % LEFT social cue condition
    idx_rightcue = find(data(:,6) == 1); % RIGHT social cue condition
    
    data_left = data(idx_leftcue, :);
    data_right = data(idx_rightcue, :);
     
    
    % Re-format data to "probit" structure for PMF fitting:
    
    
    
    
    
    % Fit PMFs to datasets:
    
    
    
    
    
    
    
    
    
    
end



















end