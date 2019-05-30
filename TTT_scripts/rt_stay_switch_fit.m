function [betas,stats] = rt_stay_switch_fit(data,transform)
% This function fits switch stay and coherence/disparity regressors to normalised reaction times in
% a 2afc rdk or cylinder choice task for single sessions of training and for combined
% sessions and tests whether resulting betas differ significantly. Data is
% analysed seperately for correct and incorrect trials,
% meaning that trial t-1 was always correct (won) and trial t was either
% correct or incorrect. The switch regressor is 1 on every trial t the subject
% made a choice in the opposite direction to trial t-1 (e.g. subject chose
% left target on trial t-1 and the right target on trial t) and 0 otherwise. The stay
% regressor is 1 on every trial t the subject made a choice in the sae
% direction as on trial t-1 (e.g. trial t-1 chose right and trial t chose
% right target). A third regressor contains the demeaned absolute (- and +
% coherences/disparities pooled) coherence/disparity of
% stimulus on trial t. These regressors are fitted to normalised RT data
% for each trial t. RT data is either normalized with log or 1/rt.

% Input:
% data = structure with size of number of sessions to be analysed
% Each session should be a matrix with rows for number of trials and
% columns 2 indicating choice left (0) or right (1), 3 being rts, 4
% coherence/disparity, and 5 indicating whether choice as correct (1) or
% incorrect (0)

% transform = either 'log' or 'div' for normalising rts

% Output
% betas = structrue with fields correct/wrong and single_sess and combined
% containg array of betas for all 3 regressors for single sessions and
% combined sessions
% stats = structure with same fields but also fields for p-values (p) and
% t-stats (t) testing whether betas for stay/switch are significantly
% different

% written by Maria Ruesseler, University of Oxford, Feb 2018



num_sessions = length(data); % get number of sessions

betas.correct.single_sess = nan(num_sessions,3); % matrix in which we save
% betas for each regressor for single sessions for trial t correct
betas.wrong.single_sess = nan(num_sessions,3);% matrix in which we save betas
% for each regressor for single sessions for trial t WRONG
stats.correct.single_sess.t  = nan(num_sessions,1); % matrix in which we
% save tstats from testing of significant difference of switch/stay regressors
% on trial t correct
stats.wrong.single_sess.t = nan(num_sessions,1);% matrix in which we save tstats
% from testing of significant difference of switch/stay regressors on trial t WRONG
stats.correct.single_sess.p  = nan(num_sessions,1);% matrix in which we save
% p-vals from testing of significant difference of switch/stay regressors on
% trial t correct
stats.wrong.single_sess.p  = nan(num_sessions,1);% matrix in which we save
% p-vals from testing of significant difference of switch/stay regressors on
% trial t WRONG

% after building regressors for each sessions we have to combine them after
% to one matrix to analyse their fit to the combined data
X_correct_all = []; % matrix in which we save regressors for trials t correct
X_wrong_all = []; % matrix for regressors for trials t wrong

% we need to do the same for the rts!
Y_correct_all = []; % matrix in which we save rts for trials t correct
Y_wrong_all = []; % matrix in which we save rts for trials t wrong


for i = 1:num_sessions % looping through sessions
    
    % get rts
    
    % normalise rts  = check distributions to make sure correct normal.
    % function is used either log or 1/rt should work
    
    switch transform % transform rts
        case 'log'
            rt_norm = log(data(i).data(:,3));
        case 'div'
            rt_norm = 1./(data(i).data(:,3));
            
    end
    
    % calculare standard deviation to excluse very high rts
    sd = std(rt_norm);
    
    % only keep trials with rts smaller than 2.5 * the standard deviation
    idx_rt_keep = rt_norm < (2.5 * sd);
    
    
    % trials we will actually analyse
    data_keep(i).data = data(i).data(idx_rt_keep,:);
    
    % get number of trials in sessions
    num_trials(i) = length(data_keep(i).data(:,1));
    
    % calculate mean coherences/disparities for demean later on
    coherences = unique(data_keep(i).data(:,4));
    mean_coherences = mean(coherences(coherences >= 0)); % because we only take absolute coherences (e.g. pool + and - coherences)
    
    
    
    % logical vector to find t-1 won left and won right vectors
    winright =  [0; (data_keep(i).data(1:num_trials(i)-1,5)== 1 & data_keep(i).data(1:num_trials(i)-1,2) == 1)]; % previous trial won right
    winleft  =  [0; (data_keep(i).data(1:end-1,5)== 1 & data_keep(i).data(1:end-1,2) == 0)]; % previous trial won left%
    
    %% stay
    %idx for stay trial t-1 correct and trial t correct for t-1 choice left and t-1
    %choice right separetely
    idx_correct_stay_right = winright == 1 & data_keep(i).data(:,5) == 1 & data_keep(i).data(:,2) == 1;
    idx_correct_stay_left = winleft == 1 & data_keep(i).data(:,5) == 1 & data_keep(i).data(:,2) == 0;
    
    
    %idx for stay trial t-1 correct and trial t WRONG for t-1 choice left and t-1
    %choice right separetely
    idx_wrong_stay_right = winright == 1 & data_keep(i).data(:,5) == 0 & data_keep(i).data(:,2) == 1;
    idx_wrong_stay_left = winleft == 1 & data_keep(i).data(:,5) == 0 & data_keep(i).data(:,2) == 0;
    %% switch
    
    % idx for switch trial t-1 correct and trial t correct for t-1 choice
    % left, trial t choice right and t-1 choice right and trial t choice
    % left separately 
    idx_correct_switch_right = winright == 1 & data_keep(i).data(:,5) == 1 & data_keep(i).data(:,2) == 0;
    idx_correct_switch_left = winleft == 1 & data_keep(i).data(:,5) == 1 & data_keep(i).data(:,2) == 1;
    
    % idx for switch trial t-1 correct and trial t WRONG for t-1 choice
    % left, trial t choice right and t-1 choice right and trial t choice
    % left separately 
    idx_wrong_switch_right = winright == 1 & data_keep(i).data(:,5) == 0 & data_keep(i).data(:,2) == 0;
    idx_wrong_switch_left = winleft == 1 & data_keep(i).data(:,5) == 0 & data_keep(i).data(:,2) == 1;
    
    %% get switch stay trials only
    
    % get indices for all trials t correct on switch and stay trials 
    idx_correct_trials = logical(  idx_correct_stay_right   + idx_correct_stay_left +...
        idx_correct_switch_right + idx_correct_switch_left);
    
    % get indices for all trials t WRONG on switch and stay trials 
    idx_wrong_trials = logical(idx_wrong_stay_right     + idx_wrong_stay_left + ...
        idx_wrong_switch_right   + idx_wrong_switch_left);
    
    % index into data and create new data matrices for all trials t correct
    % 
    switchstay_data_correct = data_keep(i).data(idx_correct_trials,:);
    
    
    % and all trials t WRONG
    switchstay_data_wrong = data_keep(i).data(idx_wrong_trials,:);
    
    %% final regressors
    
    % get all trials that 
    stay_correct = idx_correct_stay_right(idx_correct_trials) + idx_correct_stay_left(idx_correct_trials);
    
    stay_wrong = idx_wrong_stay_right(idx_wrong_trials) + idx_wrong_stay_left(idx_wrong_trials);
    
    switch_correct = idx_correct_switch_right(idx_correct_trials) + idx_correct_switch_left(idx_correct_trials);
    
    switch_wrong = idx_wrong_switch_right(idx_wrong_trials) + idx_wrong_switch_left(idx_wrong_trials);
    
    
    
    
    % pool  coherences
    idx_coh_correct = switchstay_data_correct(:,4) < 0;
    switchstay_data_correct(idx_coh_correct,4) = switchstay_data_correct(idx_coh_correct,4) .* -1;
    
    idx_coh_wrong = switchstay_data_wrong(:,4) < 0;
    switchstay_data_wrong(idx_coh_wrong,4) = switchstay_data_wrong(idx_coh_wrong,4) .* -1;
    
    % demean coherences
    switchstay_data_correct(:,4) = switchstay_data_correct(:,4) - mean_coherences;
    switchstay_data_wrong(:,4) = switchstay_data_wrong(:,4) - mean_coherences;
    
    %% glm - get betas
    
    X_correct{i} = [switchstay_data_correct(:,4),stay_correct, switch_correct];
    [betas_cor{i},~,stats_cor{i}] = glmfit(X_correct{i},rt_norm(idx_correct_trials),'normal','link','identity','constant','off');
    
    contrastmatrix = [0 1 -1];
    [pval_cor(i),tstats_cor(i)] = linhyptest(betas_cor{i}, stats_cor{i}.covb, 0, contrastmatrix, stats_cor{i}.dfe);
    
    
    X_wrong{i} = [switchstay_data_wrong(:,4),stay_wrong, switch_wrong];
    [betas_wro{i},~,stats_wro{i}] = glmfit(X_wrong{i},rt_norm(idx_wrong_trials),'normal','link','identity','constant','off');
    [pval_wro(i),tstats_wro(i)] = linhyptest(betas_wro{i}, stats_wro{i}.covb, 0, contrastmatrix, stats_wro{i}.dfe);
    
    betas.correct.single_sess(i,:) = betas_cor{i};
    betas.wrong.single_sess(i,:) = betas_wro{i};
    stats.correct.single_sess.t(i,:) = tstats_cor(i);
    stats.wrong.single_sess.t(i,:) = tstats_wro(i);
    stats.correct.single_sess.p(i,:) = pval_cor(i);
    stats.wrong.single_sess.p(i,:) = pval_wro(i);
    
    
    X_correct_all = [X_correct_all; X_correct{i}];
    X_wrong_all = [X_wrong_all; X_wrong{i}];
    Y_correct_all = [Y_correct_all; rt_norm(idx_correct_trials)];
    Y_wrong_all = [Y_wrong_all; rt_norm(idx_wrong_trials)];
    %     keyboard;
    
end % loop through sessions

% combined sessions
[betas.correct.combined,~,stats_cor_all] = glmfit(X_correct_all,Y_correct_all,'normal','link','identity','constant','off');

contrastmatrix = [0 1 -1];
[pval_cor_all,tstats_cor_all] = linhyptest(betas.correct.combined, stats_cor_all.covb, 0, contrastmatrix, stats_cor_all.dfe);


[betas.wrong.combined,~,stats_wro_all] = glmfit(X_wrong_all,Y_wrong_all,'normal','link','identity','constant','off');
[pval_wro_all,tstats_wro_all] = linhyptest(betas.wrong.combined, stats_wro_all.covb, 0, contrastmatrix, stats_wro_all.dfe);

stats.correct.combined.p = pval_cor_all;
stats.correct.combined.t = tstats_cor_all;

stats.wrong.combined.p = pval_wro_all;
stats.wrong.combined.t = tstats_wro_all;



end