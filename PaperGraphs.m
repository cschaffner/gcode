%% main module to generate histograms that appear in the paper

% created:    March 13, 2012
% by Christian Schaffner, c.schaffner@uva.nl

%% add the /lib directory to the search path
path([pwd '/lib'],path);
cd(fileparts(mfilename('fullpath')));

%% clear workspace and read in genetic code matrices
geneticcode;

%% set parameters of what we want to do

% equif flag
equif=0;

%% Theoretical polar requirement
A=Atheoreticpolar;

%% standard weights:
scoretype = 'theoretical Polar all-1 weights';

wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;

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

fixed = [1 2 3 10 11 18 19];
permutecode_subsets;
makePapergraph;


%% FH weights:
%weights from Freeland-Hurst
scoretype = 'theoretical Polar FH weights';

wtransit1=1;
wtransver1=0.5;
wtransit2=0.5;
wtransver2=0.1;
wtransit3=1;
wtransver3=1;

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

fixed = [1 2 3 10 11 18 19];
permutecode_subsets;
makePapergraph;

