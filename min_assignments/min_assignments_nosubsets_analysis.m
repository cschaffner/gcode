%% which assignments need to be fixed in 
% * addition to our 7 such that the SGC is optimal?
% * at all?

% by Christian Schaffner, c.schaffner@uva.nl
% idea by Gunnar Klau

% created: Nov 5, 2012

%% clear workspace 
geneticcode;

%% directories
scoretype = 'min_fixed_check';
solverdir=strcat(fileparts(mfilename('fullpath')),'/../solver/');
filenameinput=strcat(solverdir,scoretype, '.input');

%% weights

% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;

%weights from Freeland-Hurst
% wtransit1=1;
% wtransver1=0.5;
% wtransit2=0.5;
% wtransver2=0.1;
% wtransit3=1;
% wtransver3=1;

% implement weights:
B1=wtransit1*Btransit1 + wtransver1*Btransver1;
B2=wtransit2*Btransit2 + wtransver2*Btransver2;
B3=wtransit3*Btransit3 + wtransver3*Btransver3;
B=B1+B2+B3;

% trim the matrices to 20 x 20 (get rid of the STOP codon row / column)
Bmatrix = B(1:20,1:20);
B=Bmatrix;

B1 = B1(1:20,1:20);
B2 = B2(1:20,1:20);
B3 = B3(1:20,1:20);

A = Atheoreticpolar;

%% load results

load goodfixings_all1weights.mat

% potentially fixed blocks
fixed = [1 2 3 10 11 18 19];

% error-robustness value of the SGC
sgc=sum(sum(A .* Bmatrix));

sizefixmore=zeros(size(goodfixmore,2),1);

for i=1:size(goodfixmore,2)
    sizefixmore(i)=size(setdiff(goodfixmore{i},fixed),2);
    
    % do some sanity checks
    i
    ff=goodfixmore{i}       
    persize=20-size(ff,2);

    % do some random checks
    sgc=sum(sum(A .* B));
    for j=1:1000
        pp=randfixperm(20,ff);
        if sum(sum(A .* B(pp,pp)))<sgc
            fprintf('problem detected!')
            sgc
            pp                
            return;
        end
    end

    [cons,N]=CreateQAP(scoretype,ff,A,B);
    % run full branch and bound QAP solver
    % can also handle linear terms
    [status, result1] = system([solverdir,'qapbb < ',filenameinput,' | tail -2']);
    % returns a new permutation as result

    % reformat this new permutation
    [cost1,pos]=textscan(result1,'%d64',1);
    pnew1=textscan(result1(pos+1:end),'%2f');
    ppnew=pnew1{1}';
    val1=(cons+double(cost1{1}))/N;

    if val1<sgc
        fprintf('problem detected with qapbb!');
        ppnew
        return;
    end

 %   end
end    

% display indices of minimal sizefixmore
find(sizefixmore == min(sizefixmore))


