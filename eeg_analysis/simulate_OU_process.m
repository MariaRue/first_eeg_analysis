function [x] = simulate_OU_process(S,lbda,o,A,c,dt)
% This function simulates an OU process given the input in S. This can be
% any stimulus stream that could be used in the continuous motion task.
% This stream wiggles with some SD around 0. The output is a vector with
% the development of the decision variable over the OU process given the
% input S. 

% Input: 
% S     - input stream = driftrate over time, which is a vector over time 
% lbda  - parameter that influences the integration kernel of the OU
%         Process, negative numbers lead to stable attractors, positive
%         numbers including 0 lead to run-away processes 
% o     - dealy between Stimulus input and onset of integration kernel of
%         the OU process, for example if we have an o of 200, then the
%         first 200 sample inputs are only based on noise and input is 0
%         before Input is first sample of stimulus 
% A     - A is the amplitude by which the stimulus stream is scaled in the
%         OU process 
% c     - scales the noise of the Wiener process 
% dt    - is the time step between samples of S - should be smaller than 1, such
%         as 0.1 or 0.001

% Output is a vector (x) that tracks the decision variable of the OU
% process over time 

% Maria Ruesseler, University of Oxford 2018 



x(1) = 0;
stim_idx = 1; % count through stimulus vector
n = 1; % counts dt steps 
dt_step = 1/dt; % scales how many dt steps are necessary per stimulus step 


while stim_idx <=  length(S) % loop through stimulus - create OU process without noise
    for t = 1:dt_step % loop through dt steps per stim idx 
        if stim_idx>o % if we have offset only start with stimulus as input, 
            % if stim_idx is bigger than the offse
            I = S(stim_idx-o);
        else % otherwise input is 0 
            I = 0;
        end % if 
        
        % update x every time depending on input and noise 
        x(n+1) =  x(n) + (x(n)*lbda + I.*A) .* dt + sqrt(dt) .* c .* randn(1);
       
        n = n+1; 
        
    end % for loop 
    stim_idx = stim_idx + 1;
end % while loop 

x = x(1:dt_step:end); % scale x down to size of stimulus input 




end % function 