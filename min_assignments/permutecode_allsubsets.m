%% computes MS0 values using the subsets approach
% all possible permutations (120^4 = 207.36 * 10^6) are checked, 
% but only those with a better
% error-robustness than the SGC are kept

% in more details: 120^3 = 1728000 permutations are checked in the inner
% loop, which is executed 120 times.

% this allows to later figure out (in checkallperms_subsets.m) what happens 
% if we fix some of the blocks
% and in particular, what is the minimum amount of assignments we need to
% fix so that the SGC stays optimal.

% requires 20x20 matrices A and B,B1,B2,B3 to be in the workspace

% created: March 3, 2011
% cleaned up: Nov 3, 2012

% by Christian Schaffner, c.schaffner@uva.nl

%% clear workspace and read in genetic code matrices
geneticcode;
A=Atheoreticpolar;
scoretype = 'th_polar_FH_weights_subsets_all_permutations';

% % all-1 weights
% wtransit=1;
% wtransver=1;
% 
% wtransit1=wtransit;
% wtransver1=wtransver;
% wtransit2=wtransit;
% wtransver2=wtransver;
% wtransit3=wtransit;
% wtransver3=wtransver;


%weights from Freeland-Hurst
wtransit1=1;
wtransver1=0.5;
wtransit2=0.5;
wtransver2=0.1;
wtransit3=1;
wtransver3=1;


% geneticcode should be called beforehand
% to create 20x20 matrices A and B,B1,B2,B3

% implement weights:
B1=wtransit1*Btransit1 + wtransver1*Btransver1;
B2=wtransit2*Btransit2 + wtransver2*Btransver2;
B3=wtransit3*Btransit3 + wtransver3*Btransver3;
B=B1+B2+B3;

% trim the matrices to 20 x 20 (get rid of the STOP codon row / column)
B = B(1:20,1:20);
B1 = B1(1:20,1:20);
B2 = B2(1:20,1:20);
B3 = B3(1:20,1:20);

%% define subsets
% allowed permutation subsets from "Early Fixation of an Optimal Genetic Code" by
% Freeland et al. 
% Mol. Biol. Evol. 17(4):511-518: 2000
subset{1}=[1 6 10 17 18];
subset{2}=[2 7 11 12 19];
subset{3}=[3 4  8 13 14];
subset{4}=[5 9 15 16 20];

%subset{1}=[1:20];

%% generate permutations

fixed=[];
nrfixed=size(fixed,2);
nrnotfixed=20-nrfixed;

 if (nrfixed>0)
     error('call another function, no blocks should be fixed here');
 end

nrsubsets=size(subset,2);
subsetsize=zeros(nrsubsets,1);

for i=1:nrsubsets
    % get rid of fixed blocks
    ssubset{i}=setdiff(subset{i},fixed);
    subsetsize(i)=size(ssubset{i},2);
    
    fsubset{i}=zeros(1,20);
    fsubset{i}(1,ssubset{i})=ssubset{i};   
end

if (prod(factorial(subsetsize(1:3)))>10^7)
  error(strcat(num2str(prod(factorial(subsetsize))),' permutations will require too much memory'));    
end

tic
% start the matrix per of all permutations with the first block
per=int8(nonzeropermute(fsubset{1}));
% then "mix in" all the others
for i=2:nrsubsets-1
  nzperm=int8(nonzeropermute(fsubset{i}));
  indices=int8(ones(size(per,1),1)*[1:size(nzperm,1)]);
  indices=reshape(indices,size(indices,1)*size(indices,2),1);
  per=repmat(per,factorial(subsetsize(i)),1); % + int8(kron(nonzeropermute(fsubset{i}),ones(size(per,1),1)));
  per=per+nzperm(indices,:);
end

% the last one is kept separate
nzperm=int8(nonzeropermute(fsubset{4}));
toc


%% compute sum of all entries of the B matrices ("norm"), 
% i.e. the nb of edges in the (sub-graph)
Bnorm=sum(sum(B));
B1norm=sum(sum(B1));
B2norm=sum(sum(B2));
B3norm=sum(sum(B3));

%% compute values of standard genetic code (identity permutation)
sgc(1)=sum(sum(A .* B))/Bnorm;
sgc(2)=sum(sum(A .* B1))/B1norm;
sgc(3)=sum(sum(A .* B2))/B2norm;
sgc(4)=sum(sum(A .* B3))/B3norm;


%% scores for all permutations (implemented with for loop )

%% only keep permutations whose values are smaller than the sgc
%% and average values

vals=zeros(size(per,1),1);
means=zeros(size(nzperm,1),1);
minvals=zeros(size(per,1),1);
minper=zeros(size(per,1),20,'int8');
count=0;

% initialize timer
timing(1,0);


% loop over all permutations
for j=1:size(nzperm,1)
for i=1:size(per,1)
    p=per(i,:)+nzperm(j,:);

    % permute the B matrices
    % can be done more simply by B(p,p)

     % compute values
    vals(i)=sum(sum(A .* B(p,p)))/Bnorm;
    if ( vals(i) <= sgc(1) )
	count=count+1;
	minvals(count)=vals(i);
	minper(count,:)=p;
    end
    
    timing(1,j/size(nzperm,1));
    
end
means(j)=MEAN(vals);
end

 fname = strcat(scoretype, num2str(size(fixed,2)),'blcksfix ', num2str(size(vals(i,:),2)),'samples');
  if (suppression==1) 
      fname=strcat(fname,' suppressed');
  end
 
  save(fname,'means','sgc','count','minvals','minper');

return;

