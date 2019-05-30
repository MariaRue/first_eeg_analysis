
% Script to take in *_basic.m data array dataset and from it plot a PMF and
% correct and incorrect median RTs.
%
% Will also plot, on request, the RT distributions for a given list of
% coherence values.
%
% Created 2nd July 2016, Nela Cicmil, University of Oxford (DPAG)

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

axis([-110 110 0 1]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 18);
ylabel('Proportion rightward choices', 'FontSize', 18);
set(gca, 'FontSize', 18);

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

% NC comment (3rd July 2016): Only plotting the middle 7 coherence values
% for error trial RTs as there are not enough error trials on the largest
% coherences:
for c = 4:12
    
    incorrect_RTs = RTs(intersect(idx_incorrect, idx_coh(c).idx));
    median_RT_incorrect(1,c-3) = median(incorrect_RTs);
    lower_RT_incorrect(1,c-3) = median(incorrect_RTs) - quantile(incorrect_RTs, 0.25);
    upper_RT_incorrect(1,c-3) = quantile(incorrect_RTs, 0.75) - median(incorrect_RTs);
    
end

% Plot RT data:
figure;
hold on;
errorbar(coherences(4:12), median_RT_incorrect, lower_RT_incorrect, upper_RT_incorrect, 'o--', 'color', [0.5 0.5 0.5], 'markerfacecolor', 'w', 'LineWidth', 1.5);
errorbar(coherences, median_RT_correct, lower_RT_correct, upper_RT_correct, 'ok-', 'markerfacecolor', 'k', 'LineWidth', 1.5);

axis([-110 110 0 2000]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 18);
ylabel('Median RT with 1st & 3rd quantiles (ms)', 'FontSize', 18);
legend('Incorrect trials', 'Correct trials', 'Location', 'NorthEast');
set(gca, 'FontSize', 18);




% Plot RT correct and incorrect distributions for a given coherence value:

plot_RTdist = 1;

while plot_RTdist == 1
    plot_RTdist = input('Plot RT distribution? 1 = YES, 0 = NO ---->   ');
    
    if plot_RTdist == 1
        plot_coh = input('Type which coherence value to plot: ---->  ');
        
        coh_ID = find(coherences == plot_coh);
        
        correct_RTs = RTs(intersect(idx_correct, idx_coh(coh_ID).idx));
        incorrect_RTs = RTs(intersect(idx_incorrect, idx_coh(coh_ID).idx));
        
        % Plot the two RT distributions for correct and incorrect trials:
        figure;
        subplot(2,1,1);
        hold on;
        hist(correct_RTs, 40);
        title('Correct trials');
        xlabel('RT (ms)', 'FontSize', 16);
        ylabel('Number of trials', 'FontSize', 16);
        axis([0 2500 0 35]);
        set(gca, 'FontSize', 16);
        
        subplot(2,1,2);
        hist(incorrect_RTs, 40);
        title('Incorrect trials');
        xlabel('RT (ms)', 'FontSize', 16);
        ylabel('Number of trials', 'FontSize', 16);
        axis([0 2500 0 35]);
        
        set(gca, 'FontSize', 16);
        
    elseif plot_RTdist == 0
        fprintf('\n Bye :) \n');
        
    end
    
end








