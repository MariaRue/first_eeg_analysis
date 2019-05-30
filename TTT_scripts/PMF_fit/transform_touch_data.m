function [data] = transform_touch_data(orig_data)


%adds a 5th column to touch screen data indicating correct and false
%choices


% Maria Ruesseler, University of Oxford, 2017

% keyboard;
data = orig_data;


% choseright = data(1).data(:,2);
% right_rewarded = data(1).data(:,4)>0;
% correct = double(~xor(choseright,right_rewarded));

if isstruct(data)
    
    for i = 1: length(data)
        
        coherences = unique(data(i).data);
        vector = zeros(length(data(i).data(:,1)),1);
        data(i).data = horzcat(data(i).data,vector);
        
        bin_0 = data(i).data(:,4) == 0;
        
        for j = 1 : length(coherences)
            
            
            if coherences(j) < 0
                correct = data(i).data(:,4) == coherences(j) & data(i).data(:,2) == 0;
                data(i).data(correct,5) = 1;
                
            elseif coherences (j) == 0
                
                
                
                data(i).data(bin_0,5) = (sign(randn(sum(bin_0),1))+1) ./ 2;
                
                
            else
                
                correct = data(i).data(:,4) == coherences(j) & data(i).data(:,2) == 1;
                data(i).data(correct,5) = 1;
                
                
            end
            
        end
        
    end
    
    
    
else % combined data
    
    coherences = unique(data);
    vector = zeros(length(data(:,1)),1);
    %     data = horzcat(data,vector);
    
    bin_0 = data(:,4) == 0;
    
    for j = 1 : length(coherences)
        
        
        if coherences(j) < 0
            correct = data(:,4) == coherences(j) & data(:,2) == 0;
            data(correct,5) = 1;
            
        elseif coherences (j) == 0
            
            
            
            data(bin_0,5) = (sign(randn(sum(bin_0),1))+1) ./ 2;
            
            
        else
            
            correct = data(:,4) == coherences(j) & data(:,2) == 1;
            data(correct,5) = 1;
            
            
        end
        
    end
    
    
    
    
    
    
end







end