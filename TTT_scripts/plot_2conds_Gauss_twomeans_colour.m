function plot_fit = plot_2conds_Gauss_twomeans_colour(X, logit)

%
% This function plots the model fits predicted by the two-mean Gaussian
% model with separate means for the two baising conditions.
%
% X = the parameters; logit = the data that was used to fit the model;
%
% Nela Cicmil, 23rd December 2016, University of Oxford (DPAG)

stim_set = unique(logit.x);
g_stim = linspace(min(stim_set)-(max(stim_set/4)), max(stim_set)+(max(stim_set/4)), 100);

figure;
hold on;

% The parameters:
mu = X(1); % Mean for the right_bias condition
mu_1 = X(2); % Change in mean for the left_bias condition
sigma = X(3); % SD for both conditions


% 1: Fit to right_bias condition (expno 1 in logit dataset): 
idx_1 = find(logit.expno == 1);

for i = 1:size(idx_1, 2)
    
    data = zeros(1, logit.n(idx_1(i)) );
    data(1:logit.resps(idx_1(i))) = 1;
    std_1(i) = std(data);
    sem_1(i) = std_1(i) / sqrt(logit.n(idx_1(i))) ;
    clear data;
end

plot(stim_set, logit.resps(idx_1) ./ logit.n(idx_1), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

g_p1 = (1/2)*(1 + erf((g_stim - mu)./(sigma*sqrt(2))));
plot(g_stim, g_p1, '-', 'LineWidth', 1.75, 'Color', 'b');


% 2: Fit to left_bias condition (expno 2 in logit dataset):
idx_2 = find(logit.expno == 2);

for i = 1:size(idx_2, 2)
    
    data = zeros(1, logit.n(idx_2(i)));
    data(1:logit.resps(idx_2(i))) = 1;
    std_2(i) = std(data);
    sem_2(i) = std_2(i) / sqrt(logit.n(idx_2(i))) ;
    clear data;
end

orange = [1 0.5 0];
plot(stim_set, logit.resps(idx_2) ./ logit.n(idx_2), 'o', 'color', orange, 'MarkerFaceColor', orange);

g_p2 = (1/2)*(1 + erf((g_stim - (mu+mu_1))./((sigma)*sqrt(2))));
plot(g_stim, g_p2, '-', 'LineWidth', 1.75, 'color', orange);


% Create the output:
plot_fit.g_p1 = g_p1;
plot_fit.g_p2 = g_p2;

plot_fit.g_stim = g_stim;
plot_fit.g_prop1 = logit.resps(idx_1) ./ logit.n(idx_1);
plot_fit.g_prop2 = logit.resps(idx_2) ./ logit.n(idx_2);

plot_fit.sem_1 = sem_1;
plot_fit.sem_2 = sem_2;
plot_fit.stim_set = stim_set;


% Add SEM bars to the graph:
hold on; errorbar(plot_fit.stim_set, plot_fit.g_prop1, plot_fit.sem_1, 'o', 'Color', 'b', 'LineWidth', 1);
hold on; errorbar(plot_fit.stim_set, plot_fit.g_prop2, plot_fit.sem_2, 'o', 'Color', orange, 'LineWidth', 1);


axis([-110 110 0 1]);
xlabel('(LEFT)            % Coherence            (RIGHT)', 'FontSize', 14);
ylabel('Proportion rightward choices', 'FontSize', 14);
legend('Rightwards cue condition', '', 'Leftwards cue condition', '', 'Location', 'NorthWest');
set(gca, 'FontSize', 14);



    
end