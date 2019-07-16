function h = plotmse(ymat,dim,xsamp,varargin);

%  plotmse(ymat,dim,xsamp,varargin);
%
%  plots mean and standard error of ymat across dimension dim (default
%    1)
%
% uses nanmean, rather than mean, to ignore nans

if length(size(ymat))~=2
    error('ymat must be 2 dimensional matrix');
end

if nargin<2
    dim = 1;
end

if dim==2
    ymat = ymat';
elseif dim~=1
    error('dim must be 1 or 2');
end 

if nargin<3
    xsamp = 1:size(ymat,2);
end

mn = nanmean(ymat);
se = nanstd(ymat)./sqrt(size(ymat,1));

hold on;
h(1) = plot(xsamp,mn,'k','LineWidth',2,varargin{:});
h(2) = plot(xsamp,mn+se,'Color', [0.5, 0.5, 0.5],'LineWidth',1,varargin{:});
h(3) = plot(xsamp,mn-se,'Color', [0.5, 0.5, 0.5],'LineWidth',1,varargin{:});
hold off;