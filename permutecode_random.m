%% computation module for genetic code optimizations

% samples bign random permutations where fixed is a list of fixed blocks
% and computes a list of MS0 values where
% MS0(p) = sum_{x,y} A_{x,y} * B_{p(x),p(y)}

% requires 20x20 matrices A and B,B1,B2,B3 to be in the workspace
% also expects equif, fixed, bign to be set

% created: Feb 3, 2011
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


% if equif-flag is set, devide by block sizes
if (equif) 
    B = B ./ BlockSize(1:20,ones(1,20));
    B1 = B1 ./ BlockSize(1:20,ones(1,20));
    B2 = B2 ./ BlockSize(1:20,ones(1,20));
    B3 = B3 ./ BlockSize(1:20,ones(1,20));
    
    die('The Bst and BFH matrices should probably also be treated here');
end


%% scores for bign random permutations (implemented with for loop )
vals=zeros(4,bign);
vals_st=zeros(bign,1);
vals_FH=zeros(bign,1);

% loop over bign permutations
for i=1:bign
    
    p=randfixperm(20,fixed);
    
    % compute values
    vals(1,i)=sum(sum(A .* B(p,p)));
    vals(2,i)=sum(sum(A .* B1(p,p)));
    vals(3,i)=sum(sum(A .* B2(p,p)));
    vals(4,i)=sum(sum(A .* B3(p,p)));

    % compute values
    vals_st(i)=sum(sum(A .* Bst(p,p)));
    % compute values
    vals_FH(i)=sum(sum(A .* BFH(p,p)));
    
    % timing and status stuff:
    timing(1,i/bign);
    
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
