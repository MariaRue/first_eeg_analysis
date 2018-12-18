function [sse] = oueval(p,S,bdata,dt)

% This is a cost function for the OU process that can be used in fminsearch
% to recover the parameters of lambda (lbda), offset (o) and A (amplitude)
% of the integration kernel of an OU process. The cost function is the mean
% squared error between a model prediction (noise free) of a Stimulust (S) and the
% measured brain activity or simulated data from a noisy OU process. 

% Input: 
% p     - vector with 3 estimations, 1st = lambda, 2nd = offset and 3rd - amplitude
% S     - Stimulus input (same as the used to create bdata), could be
%         created with 'create_jumping_stimulus' 
% bdata - noise simulated or actual brain data to which we want to fit an
%         OU process to recover its parameters 
% dt    - is the time step between samples of S - should be smaller than 1, such
%         as 0.1 or 0.001 

% Output is the mean squared error estimation between model and data 

% Maria Ruesseler, University of Oxford, 2018 


lbda = -abs(p(1)); %lambda (decay if negative)

o = round(abs(p(2))); %offset
A = abs(p(3)); %amplitude
stim_idx = 1; % count through stimulus vector
dt_step = 1/dt;% scales how many dt steps are necessary per stimulus step
model(1) = 0;
n = 1; % counts dt steps

while stim_idx <=  length(S) % loop through stimulus - create OU process without noise
    for t = 1:dt_step
        if stim_idx>o
            I = S(stim_idx-o);
        else
            I = 0;
        end %i f
        
        model(n+1) =  model(n) + (model(n)*lbda + I.*A) .* dt; % OU process 
        % without noise to calculate model prediction

        n = n+1;
        
    end % for 
    stim_idx = stim_idx + 1;
end % while 


% scale x down to size of stimulus input 
model = model(1:dt_step:end-1);

 
sse = sum(abs((model - bdata))); % calculate mean squared error between model 
% prediction and actual data

end % function 