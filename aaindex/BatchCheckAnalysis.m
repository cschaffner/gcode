%% analysis of results produced by BatchCheck.m

% for how well they do with respect to MS0-value in the following settings:
% * permuting all blocks
% * fixing 6 blocks
% * fixing 6 blocks and using subsets

% analysis of results produced by BatchCheck.m

% created:    Oct 22, 2011
% by Christian Schaffner, c.schaffner@uva.nl

%% clear workspace and read in genetic code matrices
geneticcode;

%% read in results
load BatchCheckResult.mat;

% add a first column with the index number
result = [[1:size(result,1)]' result];
counter=0;
keep=[];

for i=1:size(result,1)
    if sum(result(i,2:7)<20)>4 
        counter=counter+1;
        keep=[result(i,:) ; keep];
    end
end

keep


%% set parameters of what we want to do

% nr of samples
bign = 10^5;

% equif-flag
equif=0;

% dummy weights:
wtransit1=0;
wtransver1=0;
wtransit2=0;
wtransver2=0;
wtransit3=0;
wtransver3=0;


%% loop over all aaindex1 values
looplist=keep(:,1);
looplist(1,1)=462;

%% loop over all indices in looplist
nrloops=size(looplist,1);
result=zeros(nrloops,7);

for j=1:nrloops
    ind=looplist(j);
    result(j,1)=ind;
    % skip aa values where some properties are unknown
    if sum(isnan(aaindex1(ind,:)))>0
        result(j,2:7)=NaN;
        continue;
    end
    A=mypdist(aaindex1(ind,:)') .^ 2;
    
    % all blocks    
    fixed = [];
    permutecode_random;    
    % how many codes were smaller than sgc for st weights?
    result(j,2)=sum(vals_st(:) < sgc_st);
    % how many codes were smaller than sgc for FH weights?    
    result(j,3)=sum(vals_FH(:) < sgc_FH);

    % fixing 6 blocks
    fixed = [1 2 3 10 11 18 19];
    permutecode_random;    
    % how many codes were smaller than sgc for st weights?
    result(j,4)=sum(vals_st(:) < sgc_st);
    % how many codes were smaller than sgc for FH weights?    
    result(j,5)=sum(vals_FH(:) < sgc_FH);
    
    % fixing 6 blocks and subsets
    fixed = [1 2 3 10 11 18 19];
    permutecode_subsets;    
    % how many codes were smaller than sgc for st weights?
    result(j,6)=sum(vals_st(:) < sgc_st);
    % how many codes were smaller than sgc for FH weights?    
    result(j,7)=sum(vals_FH(:) < sgc_FH);
        
    return;
    % timing and status stuff:
     if(timing(3,j/nrloops)==1)
        save('BatchCheckRefinedResult.mat','result','j');
     end

end
