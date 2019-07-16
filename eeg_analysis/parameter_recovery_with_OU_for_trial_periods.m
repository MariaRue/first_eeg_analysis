% This script is for analysing the EEG signal during trial periods (periods
% of integration). Aim is to recover the parameters of the integration
% kernel by fitting an OU process to the EEG data and trying to estimate
% lambda, amplitude and offset 

% % start fieldtrip  and add folders with .mat files with data structure from
% fieldtrip after pre-processing
addpath('/Users/maria/MATLAB-Drive/fieldtrip-master'); % fieldtrip tool box to analyse data
ft_defaults % start fieldtrip

addpath('/Users/maria/MATLAB-Drive/MATLAB/continous_rdk/analysis/eeg_analysis');


% get block sequence and number of trials for each block
addpath('/Users/maria/Documents/data/data.continous_rdk/EEG_pilot/sub001/behaviour'); %add path to folder with behavioural data and stim structures

%% load data

data = load('sub001_clean.mat');

data_final = data.data_final;

%% downsample data 
cfg = []; 
cfg.resamplefs = 100; 

data_down = ft_resampledata(cfg,data_final); 

%% load S structure from all sessions
% high freq data
S1 = load('sub001_sess017_behav.mat'); S{1} = S1.S; % disregard last trial in
% last (4th block) because right at end of session and not in EEG events - figure out how why that happened
S2 = load('sub001_sess022_behav.mat');S{2} = S2.S;
S3 = load('sub001_sess023_behav.mat');S{3} = S3.S;

%% low freq data

S1 = load('sub001_sess018_behav.mat'); S{1} = S1.S; % disregard last trial in
% last (4th block) because right at end of session and not in EEG events - figure out how why that happened
S2 = load('sub001_sess019_behav.mat');S{2} = S2.S;
S3 = load('sub001_sess020_behav.mat');S{3} = S3.S;
S4 = load('sub001_sess025_behav.mat');S{4} = S4.S;

%%

% genearate counters for saving trial indices in condition vectors for each
% block
countITIS_INTES = 0;
countITIS_INTEL = 0;
countITIL_INTES = 0;
countITIL_INTEL = 0;

% pre-allocate space for all 4 conditions
cond1 = zeros(3,2);
cond2 = zeros(3,2);
cond3 = zeros(3,2);
cond4 = zeros(3,2);
% coutner of start index of first trial for each block in each condition
trial = 1;

trial_stim = []; 
% loop through sessions
count_ch = 0;
indices = []; 
for s = 1:numel(S)
    sess = S{s}; % get session
    
    % loop through blocks of a session
    % for each block have one condition vector in which we save the trial
    % indices they correspond to in the EEG data
    blocknum = 4;
% %     for low freq 
%         if s == 4 
%             blocknum = 1;
%         else 
%             blocknum = 4;
%             
%         end
    for b = 1 : blocknum
        
        
        % get block id, 1 - 4 each number corresponding to one of the following below
        block = sess.block_ID_cells{b};
        
        % get stim indices for each trial 
        
stim = sess.coherence_frame_org{b}; % complete stimulus tream 
trials = sess.mean_coherence_org{b}; % stream with trial numbers 

indx = find(trials(1:end-1) == 0 & trials(2:end) > 0);
indx_n = find(trials(1:end-1) == 0 &trials(2:end) < 0);

indx_all = sort([indx; indx_n]); 

% get trial periods out ouf stim with 2 sec before and a couple of sec
% after 
off_b = 50;
off_a = 149; 

% %trial_periods0{s} = zeros(length(indx(1)-off_b:indx(1)+off_a ),length(indx)-1,4); 
% for i = 1:length(indx)-1
%    trial_periods = [];
% trial_periods(i,:) = stim(indx_all(i)-off_b : indx_all(i)+off_a )'; 
% count_ch = count_ch + 1; 
% %trial_periods0{s}(off_b+2:end-(off_a - 150),i,b) = 1;  
% end 

% trial_stim = [trial_stim; trial_periods]; 
%         
        % get number of trials in that block
        total_trials = max(max(sess.blocks_shuffled{b}));
        
        % switch through block ids
        switch block
            
            
            % ITIS short, INTE short
            case '1'
                %update counter for this block condition
                [indx_all] = get_stimIdx(b,sess);
                countITIS_INTES = countITIS_INTES+1;
                    if s == 3
                       total_trials = total_trials-1;
                       indx_all(end) = [];
                      
                    end
                cond1(countITIS_INTES,:) = [trial, trial+(total_trials-1)];
              indices = [indices; indx_all];
                % ITIS short, INTE long
            case '2'
                %update counter for this block condition
                 [indx_all] = get_stimIdx(b,sess);
                countITIS_INTEL = countITIS_INTEL+1;
                %
                    if s == 3
                        total_trials = total_trials -1;
                        indx_all(end) = [];
                    end
                
                cond2(countITIS_INTEL,:) = [trial, trial+(total_trials-1)];
                 indices = [indices; indx_all];
                % ITIS long, INTE short
            case '3'
                
                 [indx_all] = get_stimIdx(b,sess);
                %update counter for this block condition
                countITIL_INTES = countITIL_INTES+1;
                    if s == 1 % for high freq
                       total_trials = total_trials-1;
                       indx_all(end) = [];
                    end
                
                
%                 % for low freq
%                 if s == 2
%                     
%                     total_trials = 5;
%                     
%                 end
                
                
                cond3(countITIL_INTES,:) = [trial, trial+(total_trials-1)];
                 indices = [indices; indx_all];
                % ITIS long, INTE long
            case '4'
                 [indx_all] = get_stimIdx(b,sess);
                %update counter for this block condition
                countITIL_INTEL = countITIL_INTEL+1;
                
                cond4(countITIL_INTEL,:) = [trial, trial+(total_trials-1)];
                
                indices = [indices; indx_all];
        end % end of switch condition
        
        % update start index for next block
        trial = trial + total_trials;
        
        %trial_periods0{s} = zeros(length(indx(1)-off_b:indx(1)+off_a ),length(indx)-1,4); 
        
        length_indx = length(indx_all); 
          if s == 2 && b == 1
        
        length_indx = length_indx - 1;  
          end 
     trial_periods = [];
for i = 1:length_indx
 
  
trial_periods(i,:) = stim(indx_all(i)-off_b : indx_all(i)+off_a )'; 
count_ch = count_ch + 1; 

b
s 
%trial_periods0{s}(off_b+2:end-(off_a - 150),i,b) = 1;  
end 

trial_stim = [trial_stim; trial_periods]; 
        
    end  % end loop through blocks within a session
    
end % loop through sessions
%% recover parameters from trial periods 


% params of OU process to recover 

 dt = 0.1; 


%  define bounds of parameters lambda, offset, amplitude

bounds_l = [-1 0]; % lambda
bounds_o = [0 900]; % offset
bounds_A = [0 100]; % amplitude

% generate 10 equally spaced variables within these bounds
lambda_grid = logspace(bounds_l(1),bounds_l(2),10);
offset_grid = round(linspace(bounds_o(1),bounds_o(2),10));
amplitude_grid = linspace(bounds_A(1),bounds_A(2),10);

trial_nums = [cond1(1,1):cond1(1,2), cond1(2,1):cond1(2,2), cond1(3,1):cond1(3,2)];


            for rep = 1:length(trial_nums) % for num of trials in a condition 
                
                clear pstart 
                clear p0
                clear S 
                clear x
                sse_org = 1000000000;% initial value new calculated cost function is compared to
                % generate a stimulus and simuluate ou process
       
                rep
                
             [S]  = trial_stim(trial_nums(rep),:); 
              % S = randn(6000,1) .* 0.02; 

             
                
                % loop through all possible parameter combinations
                for l_test = 1:length(lambda_grid)
                    
                    for off_test = 1:length(offset_grid)
                        
                        for a_test = 1:length(amplitude_grid)
                            
                            p0(1) = lambda_grid(l_test);
                            p0(2)  = offset_grid(off_test);
                            p0(3) = amplitude_grid(a_test);
                            
                            sse = oueval(p0,S,data_down.trial{trial_nums(rep)}(54,:),dt);
                        
                            cost.noise.sse{off_test}(l_test,a_test,rep) = sse;
%                         if p0(1) == lambda_grid(4) && p0(2) == offset_grid(3) && p0(3) == amplitude_grid(6)
%                             keyboard; 
%                         end
                         

                            if sse < sse_org % save best estimates of p from grid search
%                                 disp(p0)
%                                 disp(sse)

                                pstart = p0;
                                sse_org = sse;
                                
                               
                            end
                            
                            
                            
                        end % loop through amplitudes
                        
                    end % loop through offsets
                    
                end % loop through lambdas
                
                % now implement fminsearch
               
                fun = @(p)oueval(p,S,data_down.trial{trial_nums(rep)}(54,:),dt); % this is the correct cost function that works
                [pnew] = fminsearch(fun,pstart);
                pnew(1) = -abs(pnew(1));
                estimates(:,rep) = pnew;
                pguess(:,rep) = pstart; 
             
            end % loop through repetitions
            
            %%  
            [h,h2] = scatterbars(estimates');

%% helper function 

function [indx_all] = get_stimIdx(b,S)

        % get block id, 1 - 4 each number corresponding to one of the following below
        block = S.block_ID_cells{b};
        
        % get stim indices for each trial 
        
stim =S.coherence_frame{b}; % complete stimulus tream 
trials = S.mean_coherence{b}; % stream with trial numbers 

indx = find(trials(1:end-1) == 0 & trials(2:end) > 0);
indx_n = find(trials(1:end-1) == 0 &trials(2:end) < 0);

indx_all = sort([indx; indx_n]); 
end % function 


%% 

% this function is from tom - plotting mean and some type of variance
function [h,h2] = scatterbars(data2plot)

jit             = (rand(size(data2plot)) - 0.5) * 0.1;
[nrows, ncolz]  = size(data2plot);
offsets         = repmat(1:ncolz,nrows,1);

figure; hold on
h2  = bar(mean(data2plot));
h   = scatter((jit(:) + offsets(:)), data2plot(:));
set(h, 'MarkerFaceColor', [0.5 0.5 0.5]);
set(h, 'MarkerEdgeColor', [1 1 1]);
end



