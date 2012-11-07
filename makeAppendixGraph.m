%% generates histograms of genetic code optimizations

% requires input scores in vals(1:4,:)
% requires sgc(1:4) scores of standard genetic code (and maybe other particular codes)
% and scoretype string

% created: Nov 3, 2012

% by Christian Schaffner, c.schaffner@uva.nl

%% set parameters

if (suppression==1)
  supp = ', STOP suppressed';
else
  supp = '';
end

% number of bins
  nbins=33;
  
  maxxx=max(max(vals(2:4,:)));
  minnn=min(min(vals(2:4,:)));
  
%% make histograms for Score0, Score1, Score2 and Score3 
  % "flexible" bins for Score0, just number is fixed
  [nn{1},xx{1}] = hist(vals(1,:), nbins); 
  
  for i=2:4  % fixed bins for Score1, Score2, Score3 in order to compare them
      [nn{i},xx{i}] = hist(vals(i,:), [minnn:(maxxx-minnn)/(1.2*nbins):maxxx]);
  end
  
  maxnn=max([nn{1} nn{2} nn{3} nn{4}]);
  
  
  % expects the bit col to be the column in which the histograms should
  % appear
  
  for i=1:4 
      subplot(4, 2, 2*i+col-1);
      hold on;
      % draw histogram:
      bar(xx{i},nn{i});
      
      % draw line for particular codes:
      plot([sgc(i), sgc(i)], [0, max(nn{i})],'--rs','LineWidth',1.5);
%       plot([gmc(i), gmc(i)], [0, max(nn{i})],'--gs');
%       plot([fhp(i), fhp(i)], [0, max(nn{i})],'--ms');
      hold off;
      
      % compute caption text:
      smaller=sum(vals(i,:) < sgc(i));
      promil=smaller/size(vals(i,:),2)*100;
      optperc=(mean(vals(i,:))-sgc(i))/(mean(vals(i,:))-min(vals(i,:))) * 100;
      % output caption text:
      if i==1
          scorecaption = ' MS_0';
      else
          scorecaption = strcat(' MScore_', num2str(mod(i-1,4)));
      end
      
      xlabel({strcat(scoretype,', ',scorecaption,supp) ; 
              strcat( num2str(smaller), ' codes (',sprintf('%.3f%%',promil) , ') <= sgc' )} );
%              strcat('sgc: ', sprintf('%2.2f',sgc(i)), ' min: ',sprintf('%.2f',min(vals(i,:))),' mean: ',sprintf('%.2f',mean(vals(i,:))),' std: ',sprintf('%.2f',std(vals(i,:))),' opt: ', sprintf('%.2f%%',optperc)) } );
      ylabel('Number of Codes');
      title(horzcat(sprintf('10^%i',log10(bign)),' samples'));
      
      % set axes:
      if (i>1)
        set(gca,'xlim',[min([sgc(2:4) gmc(2:4) fhp(2:4) min(min(vals(2:4,:)))]),max([sgc(2:4) gmc(2:4) fhp(2:4) max(max(vals(2:4,:)))])]);
        set(gca,'ylim',[0,maxnn]);
      else
        set(gca,'xlim',[min([fhp(1) gmc(1) sgc(1) vals(1,:)]),max([fhp(1) gmc(1) sgc(1) vals(1,:)])]);  
      end
  end
  

  % write graphic to output directory
  fname = strcat('PaperOutput/AppendixHistograms_', num2str(size(vals(i,:),2)),'samples');
  if (suppression==1) 
      fname=strcat(fname,' suppressed');
  end
%  set(gcf,'PaperOrientation','landscape');
  set(gcf,'PaperOrientation','portrait');
% set(gcf,'PaperPositionMode','manual');
%  set(gcf,'Position',[348 197 1000 712]);

  % write PDF
  print(fname, '-dpdf');
  % the actual figure from the paper was created by "printing the figure to
  % a ps file"
  
  % save fig
  saveas(gcf,fname);

%  save(fname,'vals','sgc','A','B1','B2','B3','B');
 