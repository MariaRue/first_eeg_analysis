function L = fit_2conds_Gauss_full(params, logit)

%
% This function aims to find the log-likelihood of the cumulative Gaussian
% fit to the input data (logit). Params is the parameter set that the model
% will use. The function fits two PMFs to the input data, which has two
% experimental conditions, right and left bias, for example from a reward
% cue or a social cue.
%
%
% The function is set up such that it can easily be used with fminsearch
% to find parameters for the maximum likelihood model fit.
%
% e.g. call the function in this way:
%
% [X, fval, exitflag,funcCount] = fminsearch(@fit_2conds_Gauss_full, params,
% options, logit)
%
%
% Nela Cicmil, 23rd December 2016, University of Oxford (DPAG)


% The parameters used in the cumulative Gaussian fits for the two
% conditions, right_bias and left_bias

mu = params(1); % Mean for the right_bias condition
mu_1 = params(2); % Change in mean for the left_bias condition

sigma = params(3); % SD for both conditions *which is assumed to be identical*


% Fit to right_bias condition (expno 1 in logit data):

idx_1 = find(logit.expno == 1);
y(idx_1) =  (1/2)*(1 + erf(([logit.x(idx_1)] - mu)./(sigma*sqrt(2))));


% Fit to left_bias condition (expno 2 in logit data):

idx_2 = find(logit.expno == 2);
y(idx_2) = (1/2)*(1 + erf(([logit.x(idx_1)] - (mu+mu_1))./(sigma*sqrt(2))));


p = y;
q = 1-y;
 

% Make sure there are no zeros (log cannot handle it)

p(find(p<1e-20)) = 1e-20;
q(find(q<1e-20)) = 1e-20;
 


% Find the negative log likelihood (NB there are many of different ways of
% doing this ... ?


% (n - resps) method
  L = -sum( (logit.resps .* log(p)) + ((logit.n - logit.resps) .* log(q)) );

  
% (1 - resps/n) method  
 %L = -sum( ((logit.resps./logit.n) .* log(p)) + ((1 - (logit.resps./logit.n)) .* log(q)) );
end








