function [mu, sd] = basic_fitpsf_cdf(probitin, varargin)
%
% The function
%       [mu, sd] = basic_fitpsf_cdf(probitin, varargin)
%
% fits a culmulative Gaussian to binomial data. The data is in a structure
% 'probit' which has:
% - probit.x = the stimulus difficulty level (for example, disparity or motion coherence) ;
% - probit.n = number of stimulus presentations of each stimulus ;
% - probit.resps = the number of responses in the +ve direction (check conventions) ;
% - probit.expno = ID number of experimental condition, for example, reward condition
%
%
% Using basic_fitpsf_cdf(probit, 'expno', n) will fit only the data for which
% probit.expno = n.
%
%
% Created 2nd July 2016, Nela Cicmil, University of Oxford (DPAG)
%
% Based on previous script "megreward_fitpsf_cdf.m".




% Make a copy of the input data.
probit = probitin;

% Dealing with the input arguments:
j = 1;
while(j < nargin)
    if(strcmpi(varargin{j},'expno'))
        % if a single expno is given, extract this data from the full set by
        % matching the expno field in the input data structure.
        j = j+1;
        expno = varargin{j};
        idx = find([probit.expno] == expno);
        probit.n = probit.n(idx);
        probit.expno = probit.expno(idx);
        probit.x = probit.x(idx);
        probit.resps = probit.resps(idx);
        j = j+1;
    end;
end;

% Calculate p: proportions of choices in +ve rotation direction.
p = [probit.resps]./[probit.n];
for j = 1:length(probit)
    probit.p(j) = p(j);
end;

% Set the fminsearch optimization options (use default)
options = optimset([]);

% If there are fewer than three points on the PSF, cannot fit the cdf.
% Exit the function with exit flag = 0.
if length(p) < 2
    results.fit(1) = NaN;
    results.fit(2) = NaN;
    results.exit = 0;
    return;
end

% Make an initial guess at the parameters to fit. x(1) is the mean, x(2) is
% the SD.

% mean:
w = exp(-(p - 0.5).^2 ./0.1);
x(1) = sum(w .* [probit.x] .* [probit.n])/sum(w .*[probit.n]);

% SD:
nsd = erfinv((p * 2) -1);

idx = find(nsd > 2);
nsd(idx) = 2;
idx = find(nsd < -2);
nsd(idx) = -2;

if(max(nsd) > 0.1)
    idx = find(abs(nsd) > 0.1);
elseif(max(nsd) > 0.02)
    idx = find(abs(nsd) > 0.02);
else
    idx = 1:length(nsd);
end

nsd(find(nsd == 0)) = mean(nsd)/10;

x(2) = sum((([probit.x(idx)] - x(1))./nsd(idx)) .* [probit.n(idx)])./sum([probit.n(idx)]);

if(x(2) > range([probit.x]))
    x(2) = range([probit.x]);
end

results.preg = polyfit([probit.x],nsd,1); % Fit a straight line to nsd as function of probit.x
x(2) = 1/results.preg(1); % ?? What is the point of replacing x(2) after spending all that time estimating it above??
a = loglike(x, probit);
b = loglike([x(1) -x(2)], probit);

if(b < a)
    x(2) = -x(2);
end

% Calclate log likelihood for the initial guess for mean and SD (held in
% array 'x').
results.initial = x;
results.initlike = loglike(x,probit);
[results.fit, results.loglike, results.exit] = fminsearch(@loglike, x, options, probit);
results.data = probit;

% Calculate the errorbars (SEM) for the data:
for i = 1:size(probit.n, 2)
    data = zeros(1, probit.n(i));
    data(1:probit.resps(i)) = 1;
    probit.std(i) = std(data);
    probit.sem(i) = probit.std(i) / sqrt(probit.n(i)) ;
    clear data;
end

% Plot the fitted curve:
stim_set = unique(probit.x);
g_stim = linspace(min(stim_set)-(max(stim_set/2)), max(stim_set)+(max(stim_set/2)), 100);

% Plot the observed data-points:

% Plot the observed data points and the fitted curve:
mu = results.fit(1); %(mean)
sd = results.fit(2); %(sd)
%   figure; subplot(2,1,1);
% % % % if expno == 1, HH = black smooth
   errorbar(stim_set, probit.resps ./ probit.n, probit.sem,'ko', 'MarkerFaceColor', 'k',  'LineWidth', 1);
   hold on;

   g_cdf = 0.5 + erf((g_stim - mu)/(sd * sqrt(2)))/2;
   plot(g_stim, g_cdf, 'k-','LineWidth', 1.5);
   hold on;

% NC comment (2nd July 2016): Uncomment the following lines to print
% n_trials at each data point: 

%    for i = 1:size(stim_set, 2)
% 
%      text(stim_set(i), probit.resps(i) ./ probit.n(i), num2str(probit.n(i)), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold'  );
%    end

% % % % elseif expno == 2
% % % %   gray = [0.4, 0.4, 0.4];
% % % %   % LL = gray dashed
% % % %   errorbar(stim_set, probit.resps ./ probit.n, probit.sem, 'o', 'color', gray, 'MarkerFaceColor', gray,  'LineWidth', 1);
% % % %   % errorbar(stim_set, probit.resps ./ probit.n, probit.sem, 'ko', 'MarkerFaceColor', 'w',  'LineWidth', 2);
% % % %   hold on;
% % % %
% % % %   g_cdf = 0.5 + erf((g_stim - mu)/(sd * sqrt(2)))/2;
% % % %   plot(g_stim, g_cdf, '-', 'color', gray, 'LineWidth', 1.5);
% % % %   % plot(g_stim, g_cdf, 'k-', 'LineWidth', 2);
% % % %   hold on;
% % % %
% % % %   blck = [0, 0, 0];
% % % %   for i = 1:size(stim_set, 2)
% % % %     text(stim_set(i)-0.0005, probit.resps(i) ./ probit.n(i), num2str(probit.n(i)), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom', 'color', gray , 'FontSize', 12, 'FontWeight', 'bold');
% % % %   end
% % % % elseif expno == 3
% % % %   % HL: = orange
% % % %   orange = [1 0.5 0];
% % % %   errorbar(stim_set, probit.resps ./ probit.n, probit.sem ./ 2, 'o', 'color', orange, 'MarkerFaceColor', orange,  'LineWidth', 1);
% % % %   hold on;
% % % %
% % % %   g_cdf = 0.5 + erf((g_stim - mu)/(sd * sqrt(2)))/2;
% % % %   plot(g_stim, g_cdf, '-', 'color', orange, 'LineWidth', 1.5);
% % % %   hold on;
% % % %
% % % %   for i = 1:size(stim_set, 2)
% % % %     text(stim_set(i)-0.0005, probit.resps(i) ./ probit.n(i), num2str(probit.n(i)), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom', 'color', orange , 'FontSize', 12, 'FontWeight', 'bold' );
% % % %   end
% % % % elseif expno == 4
% % % %   % LH = blue
% % % %   errorbar(stim_set, probit.resps ./ probit.n, probit.sem ./ 2, 'bo', 'MarkerFaceColor', 'b',  'LineWidth', 1);
% % % %   hold on;
% % % %
% % % %   g_cdf = 0.5 + erf((g_stim - mu)/(sd * sqrt(2)))/2;
% % % %   plot(g_stim, g_cdf, 'b-','LineWidth', 1.5);
% % % %   hold on;
% % % %
% % % %   for i = 1:size(stim_set, 2)
% % % %     text(stim_set(i)-0.0005, probit.resps(i) ./ probit.n(i), num2str(probit.n(i)), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom', 'color', 'b' , 'FontSize', 12, 'FontWeight', 'bold' );
% % % %   end
% % % % end

% BGC has hand-written the code for fitting the Gaussian to the
% psychometric data. This is probably because the inbuilt matlab functions
% (at least at the time) did not have the option for fitting two curves
% with the same SD. BGC has fitted the data using a likelihood-calculating
% function to evaluate the fit of the chosen parameters: 'fminsearch' is
% called to find the parameters with best-likelihood, essentially turning
% this into an MLE method for fitting. I have renamed his original 'llike'
% function to 'loglike' in my adaptation, just for clarity.
function result = loglike(params, probit)

mean = params(1);
sd = params(2);

y = ([probit.x] - mean)./sd; % first half of Gaussian norm. dist. function: y = (x - mu) / sigma.
%
% p & q: creating the norm. dist. function predicted values for all 'x' (dx) from the sd and
% means (probit) passed to this function by its caller function.
%

p = (erf(y/sqrt(2)) +1)/2;  % second half of Gaussian norm. dist. function: p = 1/2(1 + erf(y / sqrt(2) ) ).
q = 1-p;

% Can't have log(0).
% ...so BGC replaced all zeros in p & q with the minimum matlab value (a common hack, see
% Geoff Boynton's code for CSH project as well) (NC).
p(find(p<1e-20)) = 1e-20;
q(find(q<1e-20)) = 1e-20;

%
% Calculating the -ve log of the likelihood function? This would describe
% the goodness of fit of the model (NC). The  model with the largest
% likelihood value is the best-fitting model i.e. the MNE (maximum likelihood estimate).
%
% 1. Multiplying the number of responses at each dx by the log probability of
% their occurance as calculated above using the normal distribution.

% Is it valid to ignore the binomial coefficient in the equation below?
% The relevant binomial coefficient is supposed to be added to each dx binomial log-likelihood
% function calculation. The binomial coefficient for each dx is the same
% for both model fits. Since the same amount overall is added to each
% model's log-likelihood value, and these 2 model's values are substracted
% from each other to get the chi2-test variable, whether the binomial
% coefficient is included in the log-likehood calculation does not matter
% since it is only the relative difference that is used in the chi2-test.
% Hence this constant can be safely ignored for these stat purposes (NC).

result = sum( (([probit.n] - [probit.resps]) .* log(q)) + ([probit.resps] .* log(p) ) );

% 2. Making it negative:
result = -result;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitline = plotpsych(data, mean, sd, varargin)

color = 'r';
shown = 0;
j = 1;

while(j < nargin-2)
    if(strncmpi(varargin{j},'color',5))
        j = j+1;
        color = varargin{j};
    elseif strncmpi(varargin{j},'shown',4)
        shown = 1;
    end
    j = j+1;
end

if length(data) < 2
    return;
end

h = errorbar([data.x],[data.p],sqrt([data.p] .* (1 - [data.p])./[data.n]),'o');
set(h,'color',color);
set(h,'MarkerFaceColor',color,'linewidth',2);
hold on;
step = (max([data.x]) - min([data.x]))/100;
x = min([data.x]):step:max([data.x]);
y = 0.5 + erf((x - mean)/(sd * sqrt(2)))/2;
fitline = plot(x,y,'color',color,'linewidth',2);

if(shown)
    for j = 1:length(data)
        text(data(j).x+2*step,data(j).p,sprintf('%d',data(j).n),'color',color);
    end
end

