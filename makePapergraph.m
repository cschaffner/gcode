%% generates histogram of genetic code optimizations
% that appear in the paper

% requires input scores in values(1,:)
% requires sgc(1) scores of standard genetic code (and maybe other particular codes)

% created: March 13, 2012
% by Christian Schaffner, c.schaffner@uva.nl

%% set parameters

% number of bins
nbins=25;

%% make histograms for Score0

% "flexible" bins for Score0, just number is fixed
[nn{1},xx{1}] = hist(vals(1,:), nbins); 
  
clf; % clear figure 
  
hold on;
% draw histogram:
bar(xx{1},nn{1});

% set axes:
minx=min([sgc(1) vals(1,:)])-0.2;
maxx=max([sgc(1) vals(1,:)])+0.2;
set(gca,'xlim',[minx,maxx] );  

% draw line for particular codes:
plot([sgc(1), sgc(1)], [0, max(nn{1})],'--rs','LineWidth',2);
%[figx figy] = dsxy2figxy(gca, sgc(1), 0);
% annotation('textarrow',[figx figx],...
%     [0.82 0.12],'TextEdgeColor','none',...
%     'TextLineWidth',8,...
%     'FontSize',30,...
%     'String',{'SGC'},...
%     'LineWidth',4,...
%     'HeadLength',6,...
%     'HeadWidth',10,...
%     'HeadStyle','vback1',...
%     'Color',[1 0 0]);

hold off;

xlabel(xcaption);
ylabel('Number of Codes');
% title(horzcat(num2str(size(fixed,2)),' blocks fixed, ', num2str(size(vals(i,:),2)),' values'));

 
% write graphic to output directory
fname = strcat('PaperOutput/',scoretype);
set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperPositionMode','manual');
set(gcf,'Position',[348 197 1000 712]);
% write PDF
print(fname, '-dpdf');
  
% save fig
saveas(gcf,fname);

 