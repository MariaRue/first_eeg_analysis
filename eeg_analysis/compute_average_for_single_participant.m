function [data_avg] =  compute_average_for_single_participant(data,avg_flag, baseline_window, lrp_flag)
% this function computes the timelocked and baseline corrected average
% eithe across differet coherences or different conditions 

coherence = [30,40,50];

if lrp_flag
    
    if avg_flag 
        
        for i = 1:length(coherence)
        
        % find leftward trials for coherences, and rightward trials 
        idx_left = data.trialinfo(:,1) == coherence(i) + 100;
        idx_right = data.trialinfo(:,1) == coherence(i);
        
        cfg = []; 
        cfg.trials = idx_right; 
        data_coh{i,1} = ft_selectdata(cfg,data); 
        
                cfg = []; 
        cfg.trials = idx_left; 
        data_coh{i,2} = ft_selectdata(cfg,data); 
        
        
    cfg = [];
    data_avg{i,1} = ft_timelockanalysis(cfg,data_coh{i,1});
        
       
    cfg = [];
    data_avg{i,2} = ft_timelockanalysis(cfg,data_coh{i,2});
        end
        
    else 
        
        for i = 1:4
    
        % find leftward trials for coherences, and rightward trials 
        idx_left = (data.trialinfo(:,1) > 100 & data.trialinfo(:,2) == i );  
        idx_right = (data.trialinfo(:,1)  < 100 & data.trialinfo(:,2) == i);
      
        
        cfg = []; 
        cfg.trials = idx_right; 
        data_con{i,1} = ft_selectdata(cfg,data); 
        
        cfg = []; 
        cfg.trials = idx_left; 
        data_con{i,2} = ft_selectdata(cfg,data); 
        
        
        cfg = [];
        data_avg{i,1} = ft_timelockanalysis(cfg,data_con{i,1});
        
       
        cfg = [];
        data_avg{i,2} = ft_timelockanalysis(cfg,data_con{i,2});

            
        end 

    end 
    
    
else

if avg_flag 



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

end