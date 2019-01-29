% this script is for analysing some piloted eye tracking data. 
% aims: 

% 1) read in the edf file and obtain the x - y position of both eyes - do
% the make sense with regard to the calibration? 

% 2) can we align the behavioural triggers with the edf triggers?
% (upsampling of the behavioural triggers will be necessary) 

% 3) try to run the same GLM analysis we have in pilot_analysysLH on the
% x and y position of the eyes. 

% path to eyetracker 
addpath(genpath('/Users/Maria/Documents/Matlab/edf-converter')); 

addpath('/Users/Maria/Documents/data/data.continuous_rdk/eyetracker_pilot/sub003/short_session_wo_EEG'); 

%%  ---%%% read in the edf data %%%---
edf1 = Edf2Mat('s0se24.edf'); 
