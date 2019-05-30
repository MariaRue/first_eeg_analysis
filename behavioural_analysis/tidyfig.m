function tidyfig(f,fsize)

% tidyfig(figure_id,fontsize)
%
% tidies figures

if nargin<1
    f = gcf;
end
if nargin<2
    fsize = 16;
end
figure(f);
set(gca,'FontSize',fsize);
set(get(gca,'XLabel'),'FontSize',fsize);
set(get(gca,'YLabel'),'FontSize',fsize);
set(get(gca,'ZLabel'),'FontSize',fsize);
set(get(gca,'Title'),'FontSize',fsize,'FontWeight','bold');
set(gcf,'Color',[1 1 1]);
