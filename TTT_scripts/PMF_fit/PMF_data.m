function [percent_right,n,ci,rt,coherences] = PMF_data(data)

% Calculates Wilconson-s score confidence interval for rightward (ccw)
% movement (rotation), percentage of choices in that direction for each
% stim and median rts with 25 quantiles for correct and incorrect choices for each stim. Also
% returns list of coherences/disparities used and number of representation
% for each stim which is later needed for glmfit - logistic PMF fitting.
% This data is used in Plot_PMF_RT to plot PMFs and RT curves 

% Input: 
% data -  a matrix data with behavioural data with columns as specified in TTT_analysis_(eye)

% Output:
% percent_right - percentage of rightward (CCW) choices for each stim
% n - number of representations of each stim 
% confidence_interval - confidence interval for rightward (CCW) choices 
%                          for each stim calculated with Wilconson's score
% rt - matrix with median reaction times and upper and lower 25% quantiles of rts for each stim
%      column 1 = median rts for correct trials
%      column 2 = lower quantiles rts for correct trials 
%      column 3 = median rts for failed trials 
%      column 4 = lower quantile rts incorrect
%      column 5 = upper quantile rts correct 
%      colun  6 = upper quantile rts incorrect 

% coherences = list with stims 


% Maria Ruesseler, University of Oxford 2017 



coherences = unique(data(:,4)); % get list of coherences 

percent_right = zeros(size(coherences)); % vector to save percentage of rightward/CCW choices for each stim

n =  zeros(size(coherences)); % vector to save representations per stim


rt = zeros(length(coherences),6); % vector for rts 


for i = 1 : length(coherences) % loop through coherences/disparities 
    
    n(i) = sum(data(:,4)==coherences((i))); % representation per stim 
    
    percent_right(i) = (sum(data(:,4)==coherences(i) & data(:,2)== 1)./n(i)).*100; % percentage of rightward/CCW choices
    
    % calculate RTs
    keep_rt_correct = data(:,4)==coherences(i) & data(:,5) == 1; % index correct choices for given coherence/disparity
    keep_rt_false = data(:,4)==coherences(i) & data(:,5) == 0; % index failed choices for given coherence/disparity
    rt(i,1) = median(data(keep_rt_correct,3)); % median rt for correct choices 
    rt(i,3) = median(data(keep_rt_false,3)); % median rt for failed choices 
    
    % calculate lower 25 quantiles 
    rt(i,2) = rt(i,1) - quantile(data(keep_rt_correct,3),0.25);
    rt(i,4) = rt(i,3) - quantile(data(keep_rt_false,3),0.25);
    
    % calculate upper 25 quantiles 
    rt(i,5) = quantile(data(keep_rt_correct,3),0.75)-rt(i,1);
    rt(i,6) = quantile(data(keep_rt_false,3),0.75)-rt(i,3);   

end % loop through coherences/disparities 

%calculate confidence interval for percent right choices (Wilcoson's score interval) for performance for 95% confidence z1-a/2
%=1.96

z = 1.96; % z-value for 95% confidence interval
prob_correct = 0.5;  % probability of correct choice

term1 = (prob_correct.*(1-prob_correct))./n; % first term under square root
term2 = (z.^2)./(4.*(n.^2)); % second term under square root
term3 = (1+z.^2)./n; % term by which first terms get divided

ci = sqrt(term1 + term2)./term3; % vector with Wilson score confidence interval 




end % end function 
