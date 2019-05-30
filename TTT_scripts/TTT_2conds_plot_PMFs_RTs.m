function [logit_right, logit_left, mu, sd] = TTT_2conds_plot_PMFs_RTs(data, coherences, data_info)
%
% Script to take in a formatted "data" variable with two conditions (reward
% or social bias) and plot PMFs and correct median RTs for each condition.
%
% First we construct two logit structures, right_bias and left_bias, and
% send them separately to the script fitpsf_cdf_2conds.m. Then we plot the
% correct RTs separately by condition.
%
% Incorrect RTs are not plotted because there are already two correct RTs
% plotted (due to there being two conditions).
%
% Created 2nd July 2016, Nela Cicmil, University of Oxford (DPAG)
% Edited 23th December 2016.


% Separate data variable into right and left condition variables:
orig_data = data;
idx_right_cond = (data(:,5) == 1);
idx_left_cond = (data(:,5) == 2);
data_right = data(idx_right_cond, :);
data_left = data(idx_left_cond, :);


% Construct the logit structures for PMFS:
n_coh = size(coherences, 1);

% Construct the right_bias logit structure:
for c = 1:n_coh
    
    % Find all trials of the relevant coherence:
    idx_coh_right(c).idx = find(data_right(:,4) == coherences(c));
    % Put the coherence information into the logit structure:
    logit_right.n(c) = size(idx_coh_right(c).idx,1);
    logit_right.x(c) = coherences(c);
    logit_right.resps(c) = sum(data_right(idx_coh_right(c).idx,2));
    logit_right.expno(c) = 1;
    
end

% Construct the left_bias logit structure:
for c = 1:n_coh
    
    % Find all trials of the relevant coherence:
    idx_coh_left(c).idx = find(data_left(:,4) == coherences(c));
    % Put the coherence information into the logit structure:
    logit_left.n(c) = size(idx_coh_left(c).idx,1);
    logit_left.x(c) = coherences(c);
    logit_left.resps(c) = sum(data_left(idx_coh_left(c).idx,2));
    logit_left.expno(c) = 2;
    
end


% Plot the PMFs:
figure;
[mu.right, sd.right] = fitpsf_cdf_2conds(logit_right, 'expno', 1);
hold on;
[mu.left, sd.left] = fitpsf_cdf_2conds(logit_left, 'expno', 2);

axis([-110 110 0 1]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 14);
ylabel('Proportion rightward choices', 'FontSize', 14);
legend('Rightwards cue condition', '', 'Leftwards cue condition', '', 'Location', 'NorthWest');
set(gca, 'FontSize', 14);

% Write the PMF figure title:
if data_info.combined
    title(sprintf('Combined data for sessions %d to %d. Right cue: mean = %.2f; SD = %.2f; left cue: mean = %.2f; SD = %.2f.', data_info.sessions(1), data_info.sessions(2), mu.right, sd.right, mu.left, sd.left), 'FontSize', 14);
else
    title(sprintf('Data for session %d. Right cue: mean = %.2f; SD = %.2f; left cue: mean = %.2f; SD = %.2f.', data_info.session, mu.right, sd.right, mu.left, sd.left), 'FontSize', 14);
end



% Find correct trial indices, and coherence indices, for right and left bias datasets:
idx_correct_bin = zeros(1,size(data_right,1));
for a = 1:size(data_right,1)
    if data_right(a,4) == 0
        idx_correct_bin(a) = (sign(randn)+1) ./ 2;
        
    elseif ((data_right(a,2)*2 - 1) * data_right(a,4)) > 0
        idx_correct_bin(a) = 1;
        
    end
end
idx_correct_right = find(idx_correct_bin == 1);

idx_correct_bin = zeros(1,size(data_left,1));
for a = 1:size(data_left,1)
    if data_left(a,4) == 0
        idx_correct_bin(a) = (sign(randn)+1) ./ 2;
        
    elseif ((data(a,2)*2 - 1) * data_left(a,4)) > 0
        idx_correct_bin(a) = 1;
        
    end
end
idx_correct_left = find(idx_correct_bin == 1);


% Format RT datasets for right and left bias conditions:
RTs_right = data_right(:,3) * 1000; % Convert RT data to ms units.
RTs_left = data_left(:,3) * 1000; % Convert RT data to ms units.

for c = 1:n_coh % For each stimulus coherence level:
    
    correct_RTs_right = RTs_right(intersect(idx_correct_right, idx_coh_right(c).idx));
    median_RT_correct_right(1,c) = median(correct_RTs_right);
    lower_RT_correct_right(1,c) = median(correct_RTs_right) - quantile(correct_RTs_right, 0.25);
    upper_RT_correct_right(1,c) = quantile(correct_RTs_right, 0.75) - median(correct_RTs_right);
    
    correct_RTs_left = RTs_left(intersect(idx_correct_left, idx_coh_left(c).idx));
    median_RT_correct_left(1,c) = median(correct_RTs_left);
    lower_RT_correct_left(1,c) = median(correct_RTs_left) - quantile(correct_RTs_left, 0.25);
    upper_RT_correct_left(1,c) = quantile(correct_RTs_left, 0.75) - median(correct_RTs_left);
    
end


% Plot RT datasets:
orange = [1 0.5 0];
figure;
errorbar(coherences, median_RT_correct_right, lower_RT_correct_right, upper_RT_correct_right, 'ob-', 'markerfacecolor', 'b', 'LineWidth', 1.5);
hold on;
errorbar(coherences, median_RT_correct_left, lower_RT_correct_left, upper_RT_correct_left, 'o-', 'Color', orange, 'markerfacecolor', orange, 'LineWidth', 1.5);

axis([-110 110 0 2000]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 14);
ylabel('Median RT with 1st & 3rd quantiles (ms)', 'FontSize', 14);
legend('Rightwards cue condition', 'Leftwards cue condition', 'Location', 'NorthEast');
set(gca, 'FontSize', 14);

% Write the RT figure title:
if data_info.combined
    title(sprintf('Combined data for sessions %d to %d', data_info.sessions(1), data_info.sessions(2)), 'FontSize', 14);
else
    title(sprintf('Data for session %d.', data_info.session), 'FontSize', 14);
end




end








