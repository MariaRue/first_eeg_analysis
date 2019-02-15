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
path = 'eyetrack'; 
% path = 'pilot'; 
switch  path
    
    case 'pilot'
        
       


addpath('/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub004/eye'); 

filepath =  '/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub004/eye';
bhvpath = '/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub004/behaviour';

    case 'eyetrack'
       

addpath('/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub001/eyetracker_test/eyetracker'); 

filepath =  '/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub001/eyetracker_test/eyetracker';
bhvpath = '/Users/Maria/Documents/data/data.continuous_rdk/EEG_pilot/sub001/eyetracker_test/behaviour'; 
        
end 


%%  ---%%% read in the edf data and behavioural data%%%---

nsess = 2; 
subid = 0; 
session = [28,29]; 
for i  = 1:nsess 
    filename = sprintf('s%dse%d.edf',subid,session(i));
    file_to_load = fullfile(filepath,filename); 

edf{i} = Edf2Mat(file_to_load); 

bhv_file = sprintf('sub%03.0f_sess%03.0f_behav.mat',subid,session(i));
bhv_to_load = fullfile(bhvpath, bhv_file); 
bhv{i} = load(bhv_to_load);

end 

%% get events with frames for actual task 


for l = 1:nsess
     ntriggers =  length(edf{l}.Events.Messages.info); % number of triggers send during session 
   
eyetriggers{l} = zeros(length(edf{l}.Samples.time),2); % vector with nSamplesx2 with first column trigger val and second column frame
idx = 0; 
for i = 18:ntriggers % loop through triggers that are actually session related 
    Tmidx = 0;
    tcounter = 0; 
    idx = idx+1; 
    if i == ntriggers  
        keyboard; 
    end 
    
    Tridx = find(edf{l}.Events.Messages.info{i} == 'T'); % idx that divides frame and trigger number  
    while ~any(Tmidx)
    
    
    Tmidx = find(edf{l}.Samples.time == (edf{l}.Events.Messages.time(i)+tcounter));
    
    tcounter = tcounter + 1; 
    % find sample idx 
    end 
    
    
    frame = str2double(edf{l}.Events.Messages.info{i}(2:Tridx-1)); % convert frame number in double 
    
    trigger = str2double(edf{l}.Events.Messages.info{i}(Tridx+1:end)); % convert trigger in double 

    eyetriggers{l}(idx,:)  = [trigger, Tmidx]; % insert in matrix 
    
   
end 
end 
%% separate by block and annulus 

for i = 1:4
    
    annulusY{i} = []; 
    normalY{i} = []; 
    
        annulusX{i} = []; 
    normalX{i} = []; 
    
end 


for i = 1:nsess 
    keyboard; 

    % find start and end of block 
   
    eend = find(eyetriggers{i}(:,1)==210); 
    
    % for Laurences data: 
    if i == 1 
        
        
        
    end 
    
   
%     if i == 1 % for Ryans data  
%         
%         estart = [estart(1:2);2223;estart(3)];
%         end 
%            
    
    if length(estart) < 4 
        keyboard; 
    end 
    
    
    
    for block = 1:4
        
        start_idx = eyetriggers{i}(estart(block),2); 
        end_idx = eyetriggers{i}(eend(block),2); 
        
    b = str2double(bhv{i}.S.block_ID_cells{block}); 
    
    
    
    if i > 1 && mod(i,2) == 0 
        
        annulusX{b} = [annulusX{b}; edf{i}.Samples.posX(start_idx : end_idx)]; 
        annulusY{b} = [annulusY{b}; edf{i}.Samples.posY(start_idx : end_idx)]; 
         
    else 
       normalX{b} = [normalX{b}; edf{i}.Samples.posX(start_idx : end_idx)]; 
       normalY{b} = [normalY{b}; edf{i}.Samples.posY(start_idx : end_idx)];  
       
        
    end % annulus 
    
    
    
    end 
    
    
end 
%% plot data 


block_cond{1} = 'ITI short, INTE short';
block_cond{2} = 'ITI short, INTE long';
block_cond{3} = 'ITI long, INTE short';
block_cond{4} = 'ITI long, INTE long';

for i = 1:4
subplot(2,2,i)
hold on 

histogram(annulusX{i}(:),'FaceColor',[1 0 0], 'EdgeColor',[1 0 0])
histogram(normalX{i}(:),'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5)

legend('annulus on', 'annulus off')


title(block_cond{i})

hold off



end 
%%

for i = 1:4
subplot(2,2,i)
hold on 

histogram(annulusY{i}(:),'FaceColor',[1 0 0], 'EdgeColor',[1 0 0])
histogram(normalY{i}(:),'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5)

legend('annulus on', 'annulus off')


title(block_cond{i})

hold off


end 


%% get events for eye calib in Matlab with timings 


for i = 1:nsess

event_idx = 9:17; % the first 8 events are eyetracking internal and not of use for us 

idx_count = 0; 
for idx = event_idx 
    idx_count = idx_count + 1; 
    

    
    eye_calib_events{idx_count} = edf1.Events.Messages.info{idx}; % get calibration event 
    eye_calib_events_time(idx_count) = edf1.Events.Messages.time(idx); % get timing of event 
    
    eye_calib_event_idx(idx_count) = find(edf1.Samples.time == eye_calib_events_time(idx_count)); 
 
end 
end

% now plot the x and y horizontal and vertical positions for both eyes 

for p = 1:9 
  
    
    ev_t = eye_calib_event_idx(p); 
    
    figure (1)
    hold on 
    plot(edf1.Samples.posX(ev_t:ev_t+1100,1),edf1.Samples.posY(ev_t:ev_t+1100,1)); 
    hold off 
    
    figure (2) 
    hold on 
     plot(edf1.Samples.posX(ev_t:ev_t+1100,2),edf1.Samples.posY(ev_t:ev_t+1100,2)); 
    hold off 
    
end 

figure (1) 
legend('centre','leftup','leftdo','rightdo','rightup','rightmi','botmi','leftmi','topmi')
ylabel('vertical pozition of one eye') 
xlabel('horizontal position of one eye') 
title('matlab calibration check - dots appeared in every corner, centre and in the middle of all borders of screen')
tidyfig

figure (2)
legend('centre','leftup','leftdo','rightdo','rightup','rightmi','botmi','leftmi','topmi')
ylabel('vertical pozition of one eye') 
xlabel('horizontal position of one eye') 
title('matlab calibration check - dots appeared in every corner, centre and in the middle of all borders of screen')
tidyfig








%% upsample behav triggers and try to match them with eye triggers 



for block = 1:length(bhv.S.trigger_vals)
    
 
 upsampled_triggers{block} = upsample_triggers(bhv.S.trigger_vals{block},100, 1000); 
 
end 
 

%% align behavioural data with EEG data

nBlocks = length(bhv.S.trigger_vals); 
nSess = 1; 

for i = 1:nSess %loop over sessions
    eye_events = eyetriggers(find(eyetriggers(:,1))); %events in eye data
    eob = find(eye_events(:,1)==210); %end of block trigger
    
    if length(eob)~=nBlocks
        error('didn''t find 4 end of blocks');
    else
        eob = [0 eob'];
        for b = 1:nBlocks % loop over blocks
        
          % it seems like all triggers from S.triggers_val have made it
          % inot the eye file :) 
             trigger_vals_behav = bhv.S.trigger_vals{b};
            tvind = find(trigger_vals_behav); tvlist = trigger_vals_behav(tvind); % a list of all the trigger values in S.trigger vals
            
            eye_events_block = eye_events(eob(b)+1:eob(b+1)); % eye events just corresponding to this block
            tvlist_eye = eye_events_block'; % a list of all trigger values that are in the eye data
            
            if length(tvlist)<length(tvlist_eye) % are there less triggers in S.trigger_vals than in the eye list? 
                error;
            elseif tvlist(end)~=210||tvlist_eye(end)~=210 % second sanity check - is final trigger 210 in both lists? (answer is yes, i.e. 210 always found
                error;
            end
            
            %now loop backwards from end of tvlist_eeg and try to match with tvlist
            teegind = length(tvlist_eye);
            
            tind = length(tvlist);
            keep = []; %this will be a list of all entries in tvlist_eeg that we keep
            while tind>0 & teegind>0
                if tvlist(tind)==tvlist_eye(teegind)
                    keep(tind) = 1;
                    teegind = teegind - 1;
                else
                    keep(tind) = 0;
                end
                tind = tind-1;
                
              
            end % while loop  
            
            if length(tvlist(find(keep)))==length(tvlist_eye) ...
                    && all(tvlist(find(keep))==tvlist_eye') %we've matched them all!
                trigger_vals_eegmatch{b} = trigger_vals_behav;
                trigger_vals_eegmatch{b}(tvind(find(~keep))) = 0; %set the other ones to 0
            else
                fprintf('Haven''t matched session %0.0f, block %0.0f\n',i,b);
            end
            
            
        end
    end
end

%% select correct eyemovement data 
eob = []; 
eob = find(eyetriggers(:,1) == 210);   
sob = find(eyetriggers(:,1)== 11); 

for b = 1:nBlocks
    
    Eyedata{b} = [downsample(edf1.Samples.posX(sob(b):eob(b),1),10), downsample(edf1.Samples.posX(sob(b):eob(b),2),10), downsample(edf1.Samples.posY(sob(b):eob(b),1),10), downsample(edf1.Samples.posY(sob(b):eob(b),2),10)]; 
    Eyedata{b} = Eyedata{b}(1:length(bhv.S.coherence_frame{b}),:); 

end 

%% build 'sliding' GLM

tf_analysis = 0; 
clear betas
nChannels = 4;
nSess = 1;
for i = 1:nSess
    
    if i == 1
        
        nBlocks = 3;
    else
        nBlocks = 4;
    end
    
    for b = 1:nBlocks
        nLags = 150; %number of lags to test (100 lags = 1s)
        
        coherence = bhv.S.coherence_frame{b}; %vector of coherence levels for this block
        coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
        % then *all* dots move in same direction, i.e. coherence = 1
        
        
        mean_coherence = bhv.S.mean_coherence{b}; % vector of mean coherences of this block - to figure out trial periods
        
        
        %coherence = coherence(1:1000); % for piloting, delete once complete
        coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
        coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherence 'jumps'
        
        integration_start = abs([0; diff(mean_coherence)])>0; %vector of trial starts
        
        button_press = trigger_vals_eegmatch{b} == 201 |... % vector of button presses during trial periods
            trigger_vals_eegmatch{b} == 202;
        %                        trigger_vals_eegmatch{i}{b} == 205 |...
        %                        trigger_vals_eegmatch{i}{b} == 206;
        
        button_press_incoh_motion = trigger_vals_eegmatch{b} == 205 |...
            trigger_vals_eegmatch{b} == 206; % vector of button presses during intertrial periods
        
        
        trial_start = trigger_vals_eegmatch{b} == 30 |...  % get start of each trial for all coherence levels
            trigger_vals_eegmatch{b} == 40 |...
            trigger_vals_eegmatch{b} == 50 |...
            trigger_vals_eegmatch{b} == 130 |...
            trigger_vals_eegmatch{b} == 140 |...
            trigger_vals_eegmatch{b} == 150;
        
        
        % regressor for prediciton error
        coherences = [];
        coherences = bhv.S.coherence_frame{b};
        diff_coherences = diff(coherences(coherence_jump));
        diff_coherences = [coherences(1); diff_coherences]; % differnce to prev cohernce for first coherence is that coherence itself
        jump_idx = find(coherence_jump);
        coherence_level_difference = zeros(size(coherences,1),1);
        coherence_level_difference(jump_idx) = abs(diff_coherences);
        
        nF = length(coherence);
        %
        regressor_list(1).value = coherence_jump;
        regressor_list(1).nLagsBack = 100;
        regressor_list(1).nLagsForward = 150;
        regressor_list(1).name = 'coherence_jump';
        
        regressor_list(2).value = coherence_jump_level;
        regressor_list(2).nLagsBack = 100;
        regressor_list(2).nLagsForward = 150;
        regressor_list(2).name = 'coherence_jump_level';
        
        regressor_list(3).value = coherence_level_difference;
        regressor_list(3).nLagsBack = 150;
        regressor_list(3).nLagsForward = 150;
        regressor_list(3).name = 'prediction error';
        
        regressor_list(4).value = button_press;
        regressor_list(4).nLagsBack = 150;
        regressor_list(4).nLagsForward = 150;
        regressor_list(4).name = 'button_press';
        
        regressor_list(5).value = button_press_incoh_motion;
        regressor_list(5).nLagsBack = 150;
        regressor_list(5).nLagsForward = 150;
        regressor_list(5).name = 'iti button press';
        
        regressor_list(6).value = trial_start;
        regressor_list(6).nLagsBack = 50;
        regressor_list(6).nLagsForward = 500;
        regressor_list(6).name = 'trial start';
        

        
     
        Fs = 1000;
        [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs);
        
        lagged_design_matrix = lagged_design_matrix'; 
        
        for ridx = 1:length(bhv.S.coherence_frame{b}); 
            if all(~isnan(Eyedata{b}(ridx,:)))
                nanIdx(ridx) = 1; 
                
            else 
                nanIdx(ridx) = 0; 
                
            end 
        end 
        
        clean_eye_data = Eyedata{b}(logical(nanIdx),:); 
        clean_design_matrix = lagged_design_matrix(logical(nanIdx),:); 

        
   
        if tf_analysis 
         for freq = 1:size(EEGdat{i}{b},2)   
                    tmp = (pinv(lagged_design_matrix')*squeeze(EEGdat{i}{b}(:,freq,:))')';
        
        for r = 1:length(regressor_list)
            betas{freq}{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
        end
        end 
            
        else
        tmp = (pinv(clean_design_matrix)*clean_eye_data)';
        
        for r = 1:length(regressor_list)
            betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
        end
        end
    end
end



for r = 1:6
    figure;
    
   
    plotmse(squeeze(betas{r}(4,:,:)),2,time_idx(r).timebins);
  

  
  xlabel(sprintf('Influence of %s on EYE x position at time (t+X) ms',time_idx(r).name));
  
  
    tidyfig;
end




