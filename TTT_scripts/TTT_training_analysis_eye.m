function [session_data, combined_data, coherences] = TTT_training_analysis_eye(path) 

% function to read in behaviroal data from MET. Data is saved in a
% structure called session data in the field data. Each experiment date is
% assigned a row number. combined data contains the behavioral data from all 
% experiment dates and session in one matrix. Behavioural data is always
% saved in a matrix with 5 columns based on sessiondata_met.m


%Column     Information

%  1:       Trial condition - only relevant for tasks with a social or specific
%           reward condition - see TTT_training_analysis.m from Nela

%  2:       choice -  0 = leftwards choice, 1 = rightwards choice  

%  3:       Reaction time (rt)

%  4:       Coherence/Disparity 

%  5:       choice - 0 = failed choice, 1 = correct choice 


% Experiment dates directories one wants to save to a structure or matrix for
% analysis should be all saved in one folder. The path to this folder needs
% to be given as input e.g.: 

%path = '/Users/maria/Documents/data/data.nhp/TTT/behav.eye/M134.Turtle/M134.Turtle.temp/cylinder_data/';


% Maria Ruesseler, University of Oxford, November 2017 



% create empty matrix for combined session and experiment data 
combined_data = [];
coherences = [];

session_data = struct;

cd(path);

files = dir;

filenames = {files([files.isdir]==1).name}';
filenames = filenames(3:end);  % remove . and ..

% old_id enables me to combine session_data from several sessions from one
% experiment day 
old_id = 0;

for i=1:length(filenames)
    cd ([path filenames{i}]); % go into date folder 
  
    files_sess = dir; 
    filenames_sess = {files_sess([files_sess.isdir]==1).name}';
    filenames_sess = filenames_sess(3:end); % remove . and .. 
    
    for j = 1: length(filenames_sess) % looping through sesssion folders within date folders
        
        
        cd ([path,filenames{i},'/',filenames_sess{j}]);
        
        sess = load('sessdesc.mat');
        trial_num = sess.sd.trial_id; 
    
        
        session_ID(j) = sess.sd.experiment_id;

        % determine task from sessdesc tags 
        if sum(strcmp('disparity',sess.sd.tags)) > 0
            
            task = 'd'; % for disparity 
            
        elseif sum(strcmp('rdk',sess.sd.tags)) > 0
            
            task = 'c'; %for coherence
            
        end
        
        %go into trials directory 
        cd ([path,filenames{i},'/',filenames_sess{j},'/','trials']);
        
        % run sessiondata_met to extract coherence/disparity, rt, choice
        % and save to matrix 
        [data, cohlist] = sessiondata_met(trial_num,task);
        
        if session_ID(j) ~= old_id % sessions are from different experiment days
        session_data(i).data = data; 
        else % sessions are from same experiment day 
            session_data(i).data = vertcat(session_data(i).data,data);
        end
        
        
        coherences = vertcat(coherences, cohlist);
        coherences = unique(coherences);
        
    old_id = sess.sd.experiment_id;   
    end % loop through session folders 
    
    
    combined_data = vertcat(combined_data,session_data(i).data);
   
end %i=1:length(filenames) looping through date folders
    


session_ID = unique(session_ID);



end % TTT_training_analysis_eye 