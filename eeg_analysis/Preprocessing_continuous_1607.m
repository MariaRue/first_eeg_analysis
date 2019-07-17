%% specify the datasets to read in by identifying the subid
% in terminal go to the data folder containing sufolders for each
% participant --> ls > txt (make sure that list contains only direct
% subject directories)
scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
%EEGdir = '/Volumes/LaCie/data/EEG';
addpath('/Users/maria/Documents/matlab/spm12');
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
ft_defaults
reference_type = 'LM_RM';

subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28, 42];
%%
for sj = 1:length(subj_list)
        
    clearvars -EXCEPT sjj scriptdir EEGdir subj_list sj reference_type 
    
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    
    subID
    
    cd (EEGdatadir)
    
    eeg_files = dir('*set');
    
    % if raw data hasn't been transformed into eeglab set files then run
    % the following code and save raw data as set files
    if isempty({eeg_files.name})
        addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
        eeglab
        
        eeg_raw = dir('*dpa');
        nSess = length({eeg_raw.name});
        for l = 1:6
            
            
            fname_load = fullfile(EEGdatadir,...
                sprintf('sub%03.0f_sess%03.0f_eeg.cdt',subID,l));
            EEG = loadcurry(fname_load, 'CurryLocations', 'False');
            
            pop_saveset(EEG,'filename',sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,l));
        end
        
    else
        
        nSess = length({eeg_files.name});
        
    end
    
    nSess = 6; 
    
    
    % load in data from one subject or created SPM file if that hasn't
    % happened yet.
    for i = 1:nSess
        
        
        
        fname_target = fullfile(EEGdatadir,...
            sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        
        
          if exist(fname_target,'file')
             D{i} = spm_eeg_load(fname_target);
             O{i} = spm_eeg_load(fname_target);
         else
            S = [];
            
            S.dataset = fullfile(EEGdatadir,sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,i));
            
            
            S.mode = 'continuous';
            D{i} = spm_eeg_convert(S);
            
            S = [];
            S.D = D{i};
            S.fsample_new = 100;
            D{i} = spm_eeg_downsample(S);
            
            
            D{i} = spm_eeg_load(fname_target);
            
            % rereference the electrodes
            S = [];
            S.D = D{i};

            S.mode = 'write';
            
            
            S.keepothers = 1; % to keep the EOG channels
            switch reference_type
                
                
                case 'LM_RM' % rereference the 61 EEG channels to L+R mastoid average
                    
                    
                    %in our montage we have 61 EEG channels, plus the right mastoid:
                    S.montage.labelorg = D{1}.chanlabels(1:62);
                    %we keep all EEG channels, but throw out the right mastoid:
                    S.montage.labelnew = D{1}.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(62);
                    %subtract 0.5 of the right mastoid in the re-reference - see
                    %p.108 of Luck book!
                    S.montage.tra(:,62) = -0.5;
                    %get rid of the right mastoid in the new channels:
                    S.montage.tra(62,:) = [];

                    
                case 'average_reference'
                     
                    %in our montage we have 61 EEG channels. We can ignore
                    %the right mastoid, as it isn't in our equation.
                    S.montage.labelorg = D{1}.chanlabels(1:61);
                    %we keep all EEG channels
                    S.montage.labelnew = D{1}.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(61)-(1/61);
                    

                    
                otherwise 
                    error('unrecognised re-referencing type')
            end
        
            D{i} = spm_eeg_montage(S);
            
        S = [];
        S.D = D{i};
        S.band = 'bandpass';
        S.freq = [0.1 30];
        D{i} = spm_eeg_filter(S);

    end
    end
    end
    

    cd (scriptdir)
    
    %% plot different EOGs 
 
     subID = 16;
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    
    
     fname_target = fullfile(EEGdatadir,...
            sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,3));
        
        D = spm_eeg_load(fname_target); 