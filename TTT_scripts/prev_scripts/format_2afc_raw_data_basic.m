
% Function to take in a list of paths to 2afc training task datasets (*.rec
% binary format) and convert them to a 3-column array " data" : 
%   1: expno (motion coherence value)
%   2: choice (0 = left, 1 = right)
%   3: RT (ms)
%
% and an array "coherences" of the coherence values in the session. 
%
%
%
% Written 2nd July 2016, Nela Cicmil, University of Oxford (DPAG)

input(1).data_path = '~/visionlab/data/m134.turtle/20160622/m134.190.1.none.2afc_rt.mdots/m134.190.1.none.2afc_rt.mdots.rec';
input(2).data_path = '~/visionlab/data/m134.turtle/20160627/m134.192.1.none.2afc_rt.mdots/m134.192.1.none.2afc_rt.mdots.rec';
input(3).data_path = '~/visionlab/data/m134.turtle/20160628/m134.193.1.none.2afc_rt.mdots/m134.193.1.none.2afc_rt.mdots.rec';
input(4).data_path = '~/visionlab/data/m134.turtle/20160629/m134.194.1.none.2afc_rt.mdots/m134.194.1.none.2afc_rt.mdots.rec';
input(5).data_path = '~/visionlab/data/m134.turtle/20160630/m134.195.1.none.2afc_rt.mdots/m134.195.1.none.2afc_rt.mdots.rec';
input(6).data_path = '~/visionlab/data/m134.turtle/20160701/m134.196.1.none.2afc_rt.mdots/m134.196.1.none.2afc_rt.mdots.rec';

n_datasets = 6;
n_blocks = 4;

for i = 1:n_datasets

    raw_data = exp_utils_read_data_binary(input(i).data_path, 'alldata', 0);
    
    
    % Make data array for each block b:
    for b = 1:n_blocks
        
        % Find all the trial with a perceptual choice given (that is,
        % correct or incorrect choice):
        idx_trials = union(find(strcmp(raw_data.resp.type(b,:), 'correct')), find(strcmp(raw_data.resp.type(b,:), 'incorrect')));
        
        % Set up the block's data array:
        block(b).data = zeros(size(idx_trials,2), 4);
        
        % Input information into the block's data array:
        block(b).data(:,2) = ([raw_data.resp.choice{b,idx_trials}]' + 1)./2; % Choices: 0 = left, 1 = right.
        block(b).data(:,3) = [raw_data.resp.reactiontime{b,idx_trials}]'./1000; % RT (converted to seconds)
        block(b).data(:,4) = [raw_data.stim.stimulus.coherence{b,idx_trials}]'; % Coherences values

        % Convert coherence values to expno numbers in first column:
        coherences = unique([raw_data.stim.stimulus.coherence{b,idx_trials}]');
        
        for t = 1:size(idx_trials,2)
            for c = 1:size(coherences,1)
                if block(b).data(t,4) == coherences(c)
                    block(b).data(t,1) = c;
                end
            end
        end
             
   
    
    end
    
    % Concatenate data array over blocks:
    
    dataset(i).data = [block(1).data; block(2).data; block(3).data; block(4).data]; 
    clear raw_data
end

% Concatenate data array over datasets:

data = [dataset(1).data; dataset(2).data; dataset(3).data; dataset(4).data; dataset(5).data; dataset(6).data];


