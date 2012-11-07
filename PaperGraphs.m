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

%% make histogram in the main body

% Theoretical polar requirement
A=Atheoreticpolar;

% FH weights:
%weights from Freeland-Hurst
scoretype = 'theoretical_Polar_FH_weights';

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
xcaption='MS_0^{FH} value'
makePapergraph;



% standard weights:
scoretype = 'theoretical_Polar_all-1_weights';

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
figure;
xcaption='MS_0 value'
makePapergraph;

%% create histograms in the Appendix

% all-1 weights:
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

% number of samples
bign = 10^7;

% no assignments fixed
fixed = [];

% molecular distance matrix with all-1 weights
A=MatrixPeter .^ 2;
scoretype = 'molecular distance squared';

% load data if it exists, otherwise: compute and save it
if (exist('PaperOutput/moldistance_data.mat'))
    load('PaperOutput/moldistance_data.mat');
    fprintf('Loaded previous data, containing %i values\n',size(vals,2))
else
    permutecode_random;
    % save data to file for being able to just regenerate figure without
    % sampling
    save('PaperOutput/moldistance_data.mat','vals','sgc');
end

clf; % clear figure 
col = 1;
makeAppendixGraph;

% Theoretical polar requirement
A=Atheoreticpolar;
scoretype = 'updated polar requirement';

% load data if it exists, otherwise: compute and save it
if (exist('PaperOutput/th_pol_req_data.mat'))
    load('PaperOutput/th_pol_req_data.mat');
    fprintf('Loaded previous data, containing %i values\n',size(vals,2))
else
    permutecode_random;
    % save data to file for being able to just regenerate figure without
    % sampling
    save('PaperOutput/th_pol_req_data.mat','vals','sgc');
end

col = 0;
makeAppendixGraph;

return;


