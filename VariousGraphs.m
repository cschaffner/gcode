%% main module to generate histograms

% created:    Feb 7, 2011
% cleaned up: Oct 21, 2011
% by Christian Schaffner, c.schaffner@uva.nl

%% add the /lib directory to the search path
path([pwd '/lib'],path);
cd(fileparts(mfilename('fullpath')));

%% clear workspace and read in genetic code matrices
geneticcode;

%% set parameters of what we want to do

% fixed blocks
%fixed = [1 2 3 10 11 14 18 19];

%fixed = [1 2 3 10 11 18 19];
%fixed = [1 18 19];
%fixed = [1 2 10 11 18 19];
%fixed = [1 2 11 18 19];
%fixed = [1 2 18 19];

%fixed = [1 2 3 10 11 18 19 4 12 14 15];
fixed =[];

% how many samples
bign = 10^4;

% suppression of stop codons
% set in geneticcode.m

% equif flag
equif=0;

%% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;

% weights from Freeland et al.
% wtransit=1;
% wtransver=1;
% 
% wtransit1=wtransit;
% wtransver1=wtransver;
% wtransit2=wtransit;
% wtransver2=wtransver;
% wtransit3=wtransit;
% wtransver3=wtransver;

%assign weights to different changes (sub-matrices of edges)
%weights from Freeland-Hurst
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0.5;
% wtransver2=0.1;
% wtransit3=1;
% wtransver3=1;

% Gunnar's weights
% wtransit1=1;
% wtransver1=1;
% wtransit2=0;
% wtransver2=0;
% wtransit3=1;
% wtransver3=1;
% scoretype = 'Polar Gunnar weights';

% Harry's weights
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0;
% wtransver2=0;
% wtransit3=1;
% wtransver3=1;
% scoretype = 'Polar Harry weights';

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

%% Theoretical polar requirement
scoretype = 'theoretical Polar FH weights';
A=Atheoreticpolar;

CreateQAP(scoretype,fixed,A,B);
return;

permutecode_random;
makegraph;

return;

%% Original polar requirement
scoretype = 'original Polar all-1-weights';
A=Apolar;

fixed = [1 2 3 10 11 18 19];
permutecode_subsets;
makegraph;

return;


%% Wimley & White
 
scoretype = 'AWW all-1-weights';
A=AWW;


permutecode_random;
figure;
makegraph;

return;


%% Theoretical polar requirement
scoretype = 'theoretical Polar all-1-weights';
equif=0;
A=Atheoreticpolar;
permutecode_random;
makegraph;

return;


%% optimized Higgs weights
wopt=[ 0.1151         0         0         0    0.0519         0         0    0.2608    0.5721];
for i=1:9
    AA{i} = mypdist(combHiggs(:,i));
end

A=zeros(20,20);
for j=1:9
  A=A+AA{j}*wopt(j);
end

scoretype = 'optimized Higgs weights';
equif=0;
permutecode_random;
makegraph;
return;

%% my optimized weights:
wtransit1=0;
wtransver1=0;
wtransit2=0;
wtransver2=0;
wtransit3=0.9;
wtransver3=0.1;
scoretype = 'Polar optimized weights';

equif=0;
A=Apolar;
permutecode_random;
makegraph;
return;


%% Higgs's matrix
A=AHiggs;
scoretype = 'Higgs';
equif=0;

figure;
permutecode_random_split;
makegraph_split;
return;

%% Peter's matrix
A=Peter .^ 2;
scoretype = 'Peter squared';
equif=0;

figure;
permutecode_random_split;
makegraph_split;
return;


%% from Eppstein
A=ABenner;
scoretype = 'ABenner special weights';
equif=0;

w=[ 0.6050    0.3680    0.4460    0.3360    0.5820    0.5820]';
 wtransit1=w(1);
 wtransver1=w(2);
 wtransit2=w(3);
 wtransver2=w(4);
 wtransit3=w(5);
 wtransver3=w(6);

figure;
permutecode_random;
makegraph;
return;
 


%% Polar requirement with weights
% A=ASerge;
%scoretype = 'Polar';
% scoretype = 'Benner subsets';
%A=ASerge;
%

% % devide by block size (equif)
% equif=1;
% 
% % use frequencies for amino acids
% A = ASerge .* pSerge(1:20,ones(1,20));
% scoretype = 'ASerge FH weights equif aminoweights';
% 
% 
% permutecode_random_bign;
% figure;
% makegraph;
% 
% return;

%% Grantham
scoretype = 'Grantham';
A=AGrantham;
permutecode_random;
makegraph;

%% Isoelectric
scoretype = 'Isoelectric';
A=Aisoelectric;

%assign weights to different changes (sub-matrices of edges)
%weights from Freeland-Hurst
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0.5;
% wtransver2=0.1;
% wtransit3=1;
% wtransver3=1;
% scoretype = 'Isoelectric FH weights';

permutecode_random;
makegraph;

%% Hydropathy
scoretype = 'Hydropathy';
A=Ahydropathy;
permutecode_random;
makegraph;


%% second run

% fixed blocks
fixed = [1 2 3 10 11 14 18 19];
%fixed = [1 2 3 10 11 14 18 19 5 9 15 20];
%fixed = [];

%% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;


%% Polar requirement with weights
A=Apolar;
scoretype = 'Polar';

figure;
permutecode_random;
makegraph;

%% Theoretical polar requirement
scoretype = 'theoretical Polar';
A=Atheoreticpolar;
permutecode_random;
makegraph;

%% Grantham
scoretype = 'Grantham';
A=AGrantham;
permutecode_random;
makegraph;

%% Isoelectric
scoretype = 'Isoelectric';
A=Aisoelectric;

%assign weights to different changes (sub-matrices of edges)
%weights from Freeland-Hurst
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0.5;
% wtransver2=0.1;
% wtransit3=1;
% wtransver3=1;
% scoretype = 'Isoelectric FH weights';

permutecode_random;
makegraph;

%% Hydropathy
scoretype = 'Hydropathy';
A=Ahydropathy;
permutecode_random;
makegraph;


%% third run

% fixed blocks
%fixed = [1 2 3 10 11 14 18 19];
fixed = [1 2 3 10 11 14 18 19 5 9 15 20];
%fixed = [];

%% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;


%% Polar requirement with weights
A=Apolar;
scoretype = 'Polar';

figure;
permutecode_random;
makegraph;

%% Theoretical polar requirement
scoretype = 'theoretical Polar';
A=Atheoreticpolar;
permutecode_random;
makegraph;

%% Grantham
scoretype = 'Grantham';
A=AGrantham;
permutecode_random;
makegraph;

%% Isoelectric
scoretype = 'Isoelectric';
A=Aisoelectric;

%assign weights to different changes (sub-matrices of edges)
%weights from Freeland-Hurst
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0.5;
% wtransver2=0.1;
% wtransit3=1;
% wtransver3=1;
% scoretype = 'Isoelectric FH weights';

permutecode_random;
makegraph;

%% Hydropathy
scoretype = 'Hydropathy';
A=Ahydropathy;
permutecode_random;
makegraph;


