function [Prob] = winstay_even(data)
% take input data from TTT_training_analysis and calculate propabilities of

% winstay (WSt)     - previous trial correct, same choice current trial
% winshift (WSh)    - previous trial correct, different choice current trial
% looseshift (LSh)  - previous trial false, different choice current trial
% looststay (LSt)   - previous trial false, same choice current trial

% input data can be either structure with sessions or combined sessions in
% single matrix 


% Returns bar plot of probabilities and corresponding matirx 


% Maria Ruesseler, University of Oxford, 2017 

if isstruct(data) %looking at single sessions
    
    num_session = length(data); % get number of sessions 
    
    for j = 1:num_session % loop through sessions 
        
    % setting prob variables    
    WSt = 0;
    WSh = 0;
    LSh = 0;
    LSt = 0;
    
    num_trial = size(data(j).data,1); % get number of trials per session 
    
    for i = 2 :2: num_trial % loop through trials and count occurances of each condition 
        
        if data(j).data(i-1,5) == 1
            
            
            if data(j).data(i-1,2) == data(j).data(i,2)
            
            WSt = WSt +1;
            
            
            elseif data(j).data(i-1,2) ~= data(j).data(i,2)
            
            WSh = WSh +1;
            end
           
            
        elseif data(j).data(i-1,5) == 0
            
       if data(j).data(i-1,2) == data(j).data(i,2)
            
            LSt = LSt +1;
            
            
            elseif data(j).data(i-1,2) ~= data(j).data(i,2)
            
            LSh = LSh + 1;
            end
           
            
        end
        
        
    end
    
    
    % calculate probabilities 
    WSt_prob = (WSt/(WSt + WSh)).*100;
    WSh_prob = (WSh/(WSt + WSh)).*100;
    LSt_prob = (LSt/(LSt + LSh)).*100;
    LSh_prob = (LSh/(LSt + LSh)).*100;
    
    
    
     % prepare data for plotting 
     X = {
         
     'WSt'      WSt_prob
     'LSh'      LSh_prob
     'WSh'      WSh_prob
     'LSt'      LSt_prob};
     
    subplot(ceil(num_session/2),2,j);
    bar([X{:,2}]);
    set(gca,'XTickLabel',X(:,1));
    ylabel('probability %');
    title(sprintf('session %d. win-stay-loose-shift',j));
    ylim([25 75]);
     
    Prob(j).X = X;
    
    end %
    
else % combined
    
    
    WSt = 0;
    WSh = 0;
    LSh = 0;
    LSt = 0;
    
    num_trial = size(data,1);
    
    for i = 2 : num_trial
        
        if data(i-1,5) == 1 && data(i,5) == 1
            
            WSt = WSt +1;
            
        elseif data(i-1,5) == 1 && data(i,5) == 0
            
            WSh = WSh +1;
            
        elseif data(i-1,5) == 0 && data(i,5) == 0
            
            LSt = LSt + 1;
            
        elseif data(i-1,5) == 0 && data (i,5) == 1
            
            LSh = LSh + 1;
            
            
        end
        
        
    end
    
    
    
    WSt_prob = (WSt/(WSt + WSh)).*100;
    WSh_prob = (WSh/(WSt + WSh)).*100;
    LSt_prob = (LSt/(LSt + LSh)).*100;
    LSh_prob = (LSh/(LSt + LSh)).*100;
    
    
    
     
     X = {
         
     'WSt'      WSt_prob
     'LSh'      LSh_prob
     'WSh'      WSh_prob
     'LSt'      LSt_prob};
     
 Prob.X = X;
    
    bar([X{:,2}]);
    set(gca,'XTickLabel',X(:,1));
    ylabel('probability %');
    title('combined data win-stay-loose-shift');
    ylim([0 100]);

     
     
end



end