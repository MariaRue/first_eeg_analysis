function [data_avg] =  compute_average_for_single_participant(data,avg_flag, baseline_window)
% this function computes the timelocked and baseline corrected average
% eithe across differet coherences or different conditions 
if avg_flag 
coherence = [30,40,50];


for i = 1 : 3 % sort for coherences 
   
    idx_coh = data.trialinfo(:,1) == coherence(i) | data.trialinfo(:,1) == coherence(i)+100;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_coherence{i} = ft_selectdata(cfg,data); 
    
 
    

    cfg = [];
average_ERP_timelock{i} = ft_timelockanalysis(cfg,data_coherence{i});

cfg = []; 
cfg.baseline = baseline_window;
cfg.baselinetype = 'absolute';
data_avg{i} = ft_timelockbaseline(cfg,average_ERP_timelock{i});
end

else 
    
    for bl = 1 : 4 % sort for conditions
   
    idx_coh = data.trialinfo(:,2) == bl ;
    
    cfg = [];
    cfg.trials = idx_coh; 
    data_block{bl} = ft_selectdata(cfg,data); 
    
    
    

    cfg = [];


block_average_ERP_timelock{bl} = ft_timelockanalysis(cfg,data_block{bl});

cfg = []; 
cfg.baseline = baseline_window;
cfg.baselinetype = 'absolute';

data_avg{bl} = ft_timelockanalysis(cfg,block_average_ERP_timelock{bl});


end
    
    
    
    
    
end




end