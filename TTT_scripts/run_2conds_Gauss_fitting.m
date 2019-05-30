function [stat_output] = run_2conds_Gauss_fitting(logit)
    
% This function takes in a parameter set and a logit dataset (individual or
% combined), and runs the fminsearch function with the two-condition
% Gaussian fitting functions.
%
% This function compares the cumulative Gaussian fit both with and without
% an additonal parameter, mu_1, which allows biasing condition to affect
% the location of the PMF. It returns the chi-squared test pval of the
% comparison. The function also plots the result of the fit of the two-mean
% Gaussian model.
%
% Check that the starting parameters are set to something sensible.
%
% Arguments: 
% - logit: data that is evaluated for two-mean model fit
%
% 
% Nela Cicmil, 23rd December 2016, University of Oxford

% Check starting parameters:
params = [0 0 0.5]; % (1) mu, (2) mu_1, (3) SD
options = optimset([]);

% Fit the two-mean and one-mean models to the logit dataset:
[X, fval, exitflag, funcCount] = fminsearch(@fit_2conds_Gauss_full, params, options, logit);

[X_nomu1, fval_nomu1, exitflag_nomu1, funcCount_nomu1] = fminsearch(@fit_2conds_Gauss_nomu1, params([1 3]), options, logit);

plot_fit = plot_2conds_Gauss_twomeans_colour(X, logit);


% Calculating the chi-squared test, where: 
%   fval = L for full model;
%   fval_nomu1 = l for model where reward cannot affect PMF location. 
% Models differ in only one parameter, so there is 1 degree of freedom (df). 

df = 1; 
chi = 2 * (-fval + fval_nomu1);
pval = 1 - chi2cdf(chi,df); 

% Format the output:
stat_output.pval = pval;
stat_output.chi = chi; 
stat_output.X = X; 
stat_output.fval = fval; 
stat_output.exitflag = exitflag; 
stat_output.funcCount = funcCount; 
stat_output.X_nomu1 = X_nomu1; 
stat_output.fval_nomu1 = fval_nomu1; 
stat_output.exitflag_nomu1 = exitflag_nomu1;
stat_output.funcCount_nomu1 = funcCount_nomu1;


end