%% computation module for genetic code optimizations

% first compute all permutations where fixed is a list of fixed blocks
% and blocks in subsets are only permuted among them

% and computes a list of MS0 values where
% MS0(p) = sum_{x,y} A_{x,y} * B_{p(x),p(y)}

% requires 20x20 matrices A and B,B1,B2,B3 to be in the workspace
% also expects fixed to be set

% created: March 2, 2011
% cleaned: Oct 21, 2011
% by Christian Schaffner, c.schaffner@yva.nl

%% clear workspace and read in genetic code matrices

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
% allowed permutation subsets from "Early Fisaction of an Optimal Genetic Code" by
% Freeland et al. 
% Mol. Biol. Evol. 17(4):511-518: 2000
subset{1}=[1 6 10 17 18];
subset{2}=[2 7 11 12 19];
subset{3}=[3 4  8 13 14];
subset{4}=[5 9 15 16 20];

% Peter's:
% subset{1}=[1 2 6 7 10 11 12 17 18 19];
% subset{2}=[3 4  8 13 14 ];
% subset{3}=[5 9 15 16 20];

% two subsets only:
% subset{1}=[1 2 6 7 10 11 12 17 18 19];
% subset{2}=[3 4 5 8 9 13 14 15 16 20];


%subset{1}=[1:20];

%% generate permutations
% varying only over the subsets, excluding fixed blocks

nrfixed=size(fixed,2);
nrnotfixed=20-nrfixed;

% if (nrnotfixed>12)
%     error('this will require too much memory');
% end

fixed=sort(fixed);

ffixed=int8(zeros(1,20));
if (size(fixed,2)>0)
    ffixed(1,fixed)=fixed;
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

if (prod(factorial(subsetsize))>10^7)
  error(strcat(num2str(prod(factorial(subsetsize))),' permutations will require too much memory'));    
end

tic
% start the matrix per of all permutations with the fixed blocks
% and add the permutations of the first subset
per=repmat(ffixed,factorial(subsetsize(1)),1) + int8(nonzeropermute(fsubset{1}));
% then "mix in" all the others
for i=2:nrsubsets
    nzperm=int8(nonzeropermute(fsubset{i}));
    indices=int16(ones(size(per,1),1)*[1:size(nzperm,1)]);
    indices=reshape(indices,size(indices,1)*size(indices,2),1);
    per=repmat(per,factorial(subsetsize(i)),1); % + int8(kron(nonzeropermute(fsubset{i}),ones(size(per,1),1)));
    per=per+nzperm(indices,:);
end
toc


%% scores for all those permutations created above (implemented with for loop )
vals=zeros(4,size(per,1));
vals_st=zeros(size(per,1),1);
vals_FH=zeros(size(per,1),1);

% initialize timer
timing(2,0);

% loop over all permutations
for i=1:size(per,1)
    p=per(i,:);

    % compute values
    vals(1,i)=sum(sum(A .* B(p,p)));
    vals(2,i)=sum(sum(A .* B1(p,p)));
    vals(3,i)=sum(sum(A .* B2(p,p)));
    vals(4,i)=sum(sum(A .* B3(p,p)));

    % compute values
    vals_st(i)=sum(sum(A .* Bst(p,p)));
    % compute values
    vals_FH(i)=sum(sum(A .* BFH(p,p)));

    % timing stuff:
    timing(2,i/size(per,1));
    
end

%% compute sum of all entries of the B matrices ("norm"), 
% i.e. the nb of (possibly weighted) edges in the (sub-graph)
Bnorm=sum(sum(B));

% we scale all these values with the same norm, so that they stay nicely
% comparable
vals(1,:)=vals(1,:)/Bnorm;
vals(2,:)=vals(2,:)/Bnorm;
vals(3,:)=vals(3,:)/Bnorm;
vals(4,:)=vals(4,:)/Bnorm;

vals_st(:)=vals_st(:)/sum(sum(Bst));
vals_FH(:)=vals_FH(:)/sum(sum(BFH));


%% compute values of standard genetic code (identity permutation)
sgc(1)=sum(sum(A .* B))/Bnorm;
sgc(2)=sum(sum(A .* B1))/Bnorm;
sgc(3)=sum(sum(A .* B2))/Bnorm;
sgc(4)=sum(sum(A .* B3))/Bnorm;

sgc_st=sum(sum(A .* Bst))/sum(sum(Bst));
sgc_FH=sum(sum(A .* BFH))/sum(sum(BFH));

%% compute values of optimal codes

% Goldman code (pGoldman)
gmc(1)=sum(sum(A .* B(pGoldman,pGoldman)))/Bnorm;
gmc(2)=sum(sum(A .* B1(pGoldman,pGoldman)))/Bnorm;
gmc(3)=sum(sum(A .* B2(pGoldman,pGoldman)))/Bnorm;
gmc(4)=sum(sum(A .* B3(pGoldman,pGoldman)))/Bnorm;

% pFHpolar
fhp(1)=sum(sum(A .* B(pFHpolar,pFHpolar)))/Bnorm;
fhp(2)=sum(sum(A .* B1(pFHpolar,pFHpolar)))/Bnorm;
fhp(3)=sum(sum(A .* B2(pFHpolar,pFHpolar)))/Bnorm;
fhp(4)=sum(sum(A .* B3(pFHpolar,pFHpolar)))/Bnorm;

return;
