% barplot with raw datapoints
% usage: h = scatterbars(data)
% where data is an n x m matrix of data
% plots m bars, with bar height 'nanmean(data)'
% optional second argument jitterWidth controls x-axis spread of dots
%
% example usage
% x         = rand(30,4);       % generate some raw data
% x(1:20,3) = NaN;              % scatterbars can handle missing data
% figure; h = scatterbars(x);   % plot

function [h] = scatterbars(data, jitterWidth)

% specify dot jitter width if not already specified
if nargin < 2
    jitterWidth = 0.2;
end

% random jitter in the x-dimension for the dots
jitter          = (rand(size(data)) - 0.5) * jitterWidth;
[nRows, nCols]  = size(data);
offsets         = repmat(1:nCols,nRows,1);

hold on

% bar heights
means           = nanmean(data);

% plot
for i = 1:nCols
    % bars
    h.bar(i)    = bar(i, means(i));
    % dots
    h.dots(i)   = scatter((jitter(:,i) + offsets(:,i)), data(:,i));
end

% apply some global properties to all dots
set(h.dots, 'MarkerFaceColor', [0.5 0.5 0.5]);
set(h.dots, 'MarkerEdgeColor', [1 1 1]);