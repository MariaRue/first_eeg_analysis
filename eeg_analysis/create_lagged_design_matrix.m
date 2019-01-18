function [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs)

% [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs)
% 
% INPUTS:
% regressor_list is a structure with fields for each regressor:
%   .value [1 * nEEGSamples vector of regressor's value]
%   .nLagsBack [how many lags backward does this regressor go?]
%   .nLagsForward [how many lags forward does this regressor go?]
%   .name [what is this regressor's name?]
%
% Fs is sampling frequency in Hz
%
% OUTPUTS:
% lagged_design_matrix is a design matrix for the EEG data, with rows (nLagsBack+nLagsForward+1) for each regressor and columns nEEGSamples
% time_idx is a structure with fields for each regressor:
%    .timebins [what are the timebins for this regressor, in ms, for plotting etc]
%    .dm_row_idx [which rows in the design matrix is this regressor?]
%    .name


% validate regressor_list
nReg = length(regressor_list);
nSamples = length(regressor_list(1).value);
for i = 1:nReg
    if any(size(regressor_list(i).value)~=[nSamples 1])
        error('regressor_list.value must be nEEGsamples * 1 for all regressors');
    elseif round(regressor_list(i).nLagsBack)~=regressor_list(i).nLagsBack || regressor_list(i).nLagsBack<0
        error('nLagsBack must be positive integer or 0');
    elseif round(regressor_list(i).nLagsForward)~=regressor_list(i).nLagsForward || regressor_list(i).nLagsForward<0
        error('nLagsForward must be positive integer or 0');
    end
end

%now build design matrix and time_idx
row_start = 1; %indicator for where we are in the design matrix rows

for i = 1:nReg    
    time_idx(i).name = regressor_list(i).name;
    time_idx(i).timebins = (1000/Fs) * [-regressor_list(i).nLagsBack:regressor_list(i).nLagsForward]; %time index for this regressor, in milliseconds
    
    time_idx(i).dm_row_idx = row_start:row_start+regressor_list(i).nLagsBack+regressor_list(i).nLagsForward;
    row_start = time_idx(i).dm_row_idx(end) + 1;
    
    nRows = length(time_idx(i).dm_row_idx); %number of rows for this regressor
    
    reg_zp = [zeros(regressor_list(i).nLagsForward,1); regressor_list(i).value; ...
        zeros(regressor_list(i).nLagsBack,1)]; %zero pad regressor
    
    lagged_regressor = nan(nRows,nSamples);
    for j = 1:nRows
        rowidx = nRows - j + 1;
        lagged_regressor(rowidx,:) = reg_zp(j:(j-1+nSamples));
    end
    
    lagged_design_matrix(time_idx(i).dm_row_idx,:) = lagged_regressor;
end

