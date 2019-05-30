
function [mu, sd] = TTT_basic_plot_PMF_RTs(data, coherences, data_info)
% Script to take in a formatted "data" variable and plot a PMF and
% correct and incorrect median RTs.
%
%
% Created 2nd July 2016, Nela Cicmil, University of Oxford (DPAG)
% Edited 20th December 2016.

orig_data = data;
n_coh = size(coherences, 1);

% Construct the logit structure:
for c = 1:n_coh
    
    % Find all trials of the relevant coherence:
    idx_coh(c).idx = find(data(:,4) == coherences(c));
    % Put the coherence information into the logit structure:
    logit.n(c) = size(idx_coh(c).idx,1);
    logit.x(c) = coherences(c);
    logit.resps(c) = sum(data(idx_coh(c).idx,2));
    logit.expno(c) = 1;
    
end


% Plot the PMF:
figure;
[mu, sd] = basic_fitpsf_cdf(logit, 'expno', 1);

axis([-0.04 0.04 0 1]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 14);
ylabel('Proportion rightward choices', 'FontSize', 14);
set(gca, 'FontSize', 14);

% Write the PMF figure title:
if data_info.combined
    title(sprintf('Combined data for sessions %d to %d. Mean = %.2f; SD = %.2f.', data_info.sessions(1), data_info.sessions(2), mu, sd), 'FontSize', 14);
else
    title(sprintf('Data for session %d. Mean = %.2f; SD = %.2f.', data_info.session, mu, sd), 'FontSize', 14);
end

% Find correct and incorrect trial indices, and coherence trial indices:
idx_correct_bin = zeros(1,size(data,1));

for a = 1:size(data,1)
    if data(a,4) == 0
        idx_correct_bin(a) = (sign(randn)+1) ./ 2;
        
    elseif ((data(a,2)*2 - 1) * data(a,4)) > 0
        idx_correct_bin(a) = 1;
        
    end
end
idx_correct = find(idx_correct_bin == 1);
idx_incorrect = find(idx_correct_bin == 0);


% Format RT data:
RTs = data(:,3) * 1000; % Convert RT data to ms units.

for c = 1:n_coh % For each stimulus coherence level:
    
    correct_RTs = RTs(intersect(idx_correct, idx_coh(c).idx));
    median_RT_correct(1,c) = median(correct_RTs);
    lower_RT_correct(1,c) = median(correct_RTs) - quantile(correct_RTs, 0.25);
    upper_RT_correct(1,c) = quantile(correct_RTs, 0.75) - median(correct_RTs);
    
end

% NC comment (21st December 2016): Excluding the n_exc largest stimulus
% values at each end of the range for incorrect RTs, because there are so
% few trials and the values are too noisy.

n_exc = 3; % Exclude the largest coherence values at each end of the range.

inc_start_idx = n_exc + 1;
inc_end_idx = n_coh - n_exc;

for c = inc_start_idx:inc_end_idx
    
    incorrect_RTs = RTs(intersect(idx_incorrect, idx_coh(c).idx));
    median_RT_incorrect(1,c-n_exc) = median(incorrect_RTs);
    lower_RT_incorrect(1,c-n_exc) = median(incorrect_RTs) - quantile(incorrect_RTs, 0.25);
    upper_RT_incorrect(1,c-n_exc) = quantile(incorrect_RTs, 0.75) - median(incorrect_RTs);
    
end

% Plot RT data:
figure;
hold on;
errorbar(coherences(inc_start_idx:inc_end_idx), median_RT_incorrect, lower_RT_incorrect, upper_RT_incorrect, 'o--', 'color', [0.5 0.5 0.5], 'markerfacecolor', 'w', 'LineWidth', 1.5);
errorbar(coherences, median_RT_correct, lower_RT_correct, upper_RT_correct, 'ok-', 'markerfacecolor', 'k', 'LineWidth', 1.5);

%axis([-110 110 0 2000]);
axis([-0.04 0.04 150 300]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 14);
ylabel('Median RT with 1st & 3rd quantiles (ms)', 'FontSize', 14);
legend('Incorrect trials', 'Correct trials', 'Location', 'NorthEast');
set(gca, 'FontSize', 14);

% Write the RT figure title:
if data_info.combined
    title(sprintf('Combined data for sessions %d to %d', data_info.sessions(1), data_info.sessions(2)), 'FontSize', 14);
else
    title(sprintf('Data for session %d.', data_info.session), 'FontSize', 14);
end




end








