%% latex-output of results produced by BatchCheck.m

% for how well they do with respect to MS0-value in the following settings:
% * permuting all blocks
% * fixing 7 blocks
% * fixing 7 blocks and using subsets

% analysis of results produced by BatchCheck.m

% created:    Nov 18, 2012
% by Christian Schaffner, c.schaffner@uva.nl

%% clear workspace and read in genetic code matrices
geneticcode;

%% read in results
load BatchCheckRefinedResult10_6.mat

%% produce latex output

sortresult = sortrows(result,[2 5 7]);
ranks=zeros(size(sortresult));

% compute ranks
for j=2:size(sortresult,2)
    ranks(:,j)=StandardCompetitionRankings(sortresult(:,j));
end

fid = fopen('batchcheck_table.tex','w');
for i=1:size(sortresult,1)
    for j=2:size(sortresult,2)
        fprintf(fid,'%i & (%i) &',sortresult(i,j),ranks(i,j));
    end
    fprintf(fid,'%s\\\\\n',aaind_desc{sortresult(i,1)});
end
fclose(fid);


