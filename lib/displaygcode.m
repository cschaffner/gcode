% nicely displays the according genetic code
% takes a set of 20 amino acid values, a normalization flag 
% and (optionally) a permutation

% takes as input a set of 20 values,
% a flag normalize,
% a row vector describing a permutation (if not specified,
% the identity permutation, i.e. the standard genetic code, is used)

% if the flag normalize is non-zero, the values are normalized to have
% mean 0 and standard deviation 10.

% output: displays the according genetic code

% created: 25 March 2011
% cleaned: 27 October 2011
% by Christian Schaffner, c.schaffner@uva.nl


function displaygcode(values, normalize, varargin)

%% define genetic code
  % names of amino acids
  aminos = {'Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly', 'STOP'};

  % Block structure of genetic code
    % U=0, C=1, A=2, G=3

    Block{1}= [0 0 0 ; 0 0 1];
    Block{2}= [0 0 2 ; 0 0 3 ; 1 0 0 ; 1 0 1 ; 1 0 2 ; 1 0 3];
    Block{3}= [2 0 0 ; 2 0 1 ; 2 0 2 ];
    Block{4}= [2 0 3 ];
    Block{5}= [3 0 0 ; 3 0 1 ; 3 0 2 ; 3 0 3];
    Block{6}= [0 1 0 ; 0 1 1 ; 0 1 2 ; 0 1 3 ; 2 3 0 ; 2 3 1];
    Block{7}= [1 1 0 ; 1 1 1 ; 1 1 2 ; 1 1 3];
    Block{8}= [2 1 0 ; 2 1 1 ; 2 1 2 ; 2 1 3];
    Block{9}= [3 1 0 ; 3 1 1 ; 3 1 2 ; 3 1 3];
    Block{10}=[0 2 0 ; 0 2 1];
    Block{11}=[1 2 0 ; 1 2 1];
    Block{12}=[1 2 2 ; 1 2 3];
    Block{13}=[2 2 0 ; 2 2 1];
    Block{14}=[2 2 2 ; 2 2 3];
    Block{15}=[3 2 0 ; 3 2 1];
    Block{16}=[3 2 2 ; 3 2 3];
    Block{17}=[0 3 0 ; 0 3 1];
    Block{18}=[0 3 3 ];
    Block{19}=[1 3 0 ; 1 3 1 ; 1 3 2 ; 1 3 3 ; 2 3 2 ; 2 3 3];
    Block{20}=[3 3 0 ; 3 3 1 ; 3 3 2 ; 3 3 3];
    Block{21}=[0 2 2 ; 0 2 3 ; 0 3 2 ];  % STOP codon

    % fill in code and compute sizes of blocks
    Code=zeros(4,4,4);
    BlockSize=zeros(21,1);
    for r=1:21
        BlockSize(r)=size(Block{r},1);
        for i=1:BlockSize(r)
            Code(Block{r}(i,1)+1,Block{r}(i,2)+1,Block{r}(i,3)+1)=r;
        end
    end

%% process permutation input
  optargin = size(varargin,2);      
  if (optargin>1)
      error('function should be called with a row vector of values, and possibly a row vector of permutation');
  elseif (optargin==1)
      p=varargin{1};
  elseif (optargin==0)
      p=1:21;
  end
  
  if (size(p,1)>1)
      error('argument needs to be a row vector');
  elseif (size(p,2)<20 || size(p,2)>21 )
      error('row vector needs to contain 20 or 21 distinct elements');
  end

  if (size(p,2)==20)
      % add the identity on the stop codon
      p = [p 21];
  end
  
  if (sum(p)~=210+21)
      error('elements are not distinct');
  end
  
  
%% process value input
  if normalize
      values=10*(values-mean(values))/std(values);
      formatstring='%4s %+6.1f \t'; % use fixed-point in this case
  else
      formatstring='%4s %+7.6g \t'; % use relevant digits here
  end
  aa_values=[values NaN];

%% preprocessing and output

  %reshape code to usual display
  reCode=reshape(reshape(Code,16,4)',16,4);
  
  % output
  for i=1:16
      for j=1:4
        fprintf(formatstring,aminos{p(reCode(i,j))},aa_values(p(reCode(i,j))));       
      end    
      fprintf('\n');
    if (mod(i,4)==0)
        fprintf('\n');
    end
  end
 
end
