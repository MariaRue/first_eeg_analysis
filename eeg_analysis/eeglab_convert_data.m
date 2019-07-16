% EEGdatadir = '/Users/maria/Documents/data/data.continuous_rdk/EEG_pilot/sub000/eeg/';
addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
eeglab
EEGdatadir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG','sub024','EEG');

subID = 24; 
nSess = 6;
cd(EEGdatadir)

for i = 1:nSess
clear EEG; 
    
    fname_load = fullfile(EEGdatadir,...
        sprintf('sub%03.0f_sess%03.0f_eeg.cdt',subID,i));
EEG = loadcurry(fname_load, 'CurryLocations', 'False');

pop_saveset(EEG,'filename',sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,i));
end 