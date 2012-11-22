%% which assignments need to be fixed in 
% * addition to our 7 such that the SGC is optimal?
% * at all?

% by Christian Schaffner, c.schaffner@uva.nl
% idea by Gunnar Klau

% created: Nov 5, 2012

%% clear workspace 
geneticcode;

%% weights

% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;

%weights from Freeland-Hurst
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
Bmatrix = B(1:20,1:20);

B1 = B1(1:20,1:20);
B2 = B2(1:20,1:20);
B3 = B3(1:20,1:20);

equif=0;
% if equif-flag is set, devide by block sizes
if (equif) 
    B = B ./ BlockSize(1:20,ones(1,20));
    B1 = B1 ./ BlockSize(1:20,ones(1,20));
    B2 = B2 ./ BlockSize(1:20,ones(1,20));
    B3 = B3 ./ BlockSize(1:20,ones(1,20));
    
    die('The Bst and BFH matrices should probably also be treated here');
end

scoretype = 'min_assignments';
solverdir=strcat(fileparts(mfilename('fullpath')),'/../solver/');
filenameinput=strcat(solverdir,scoretype, '.input');

%% test

%fixed = [1 2 3 10 11 18 19];
fixed = [];
nonfix=setdiff(1:20,fixed);

A = Atheoreticpolar;

goodcount=0;
bbcount=0;
skipset=0;
goodfixmore={};

% error-robustness value of the SGC
sgc=sum(sum(A .* Bmatrix));

for k=1:18
    % create a matrix where in each row contains a different permutation with
    % k out of the nonfix positions
    fixmore=mycombnk(nonfix,k);
    persize=20-size(fixed,2)-k;


    timing(4,0);
    for i=1:size(fixmore,1)        
        timing(4,i/size(fixmore,1));
        
        fprintf('now testing to fix: (already found: %i, bb runs: %i)',goodcount,bbcount);
        fixmore(i,:)
        
        ff=sort([fixed,fixmore(i,:)]);
        % check if the set ff is a superset of one we know that works
        for j=1:goodcount
            intersec=intersect(ff,goodfixmore{j});
            % if the size of the intersection is as large as the size of
            % goodfixmore{j}, then ff must be a superset of goodfixmore{j}
            if size(intersec,2)==size(goodfixmore{j},2)
                % set a flag to skip that set
                skipset=1;
                break
            end
        end
        % skip this set and continue with the next
        if skipset==1
            % reset the flag
            skipset=0;
            continue;
        end
        
        [cons,N]=CreateQAP(scoretype,ff,A,Bmatrix);
        % run heuristics
        % "qapsim.f" from http://www.seas.upenn.edu/qaplib/codes.html
        % this is the faster one, the other heuristic will only be run if
        % necessary
        [status, result1] = system([solverdir,'mysimqap < ',filenameinput,' | tail -2']);
        % returns a new permutation as result

        % "gqapd.f" GRASP algorithm
        % Notice that the other heuristic cannot handle linear terms!!!
        
        % reformat this new permutation
        [cost1,pos]=textscan(result1,'%d64',1);
        pnew1=textscan(result1(pos+1:end),'%2f');
        ppnew=pnew1{1}';
        val1=(cons+double(cost1{1}))/N;

        % we compare the obtained error-robustness with the sgc value
        % instead of looking at the actual permutation, because there are 
        % non-identical permutations which give the sgc as well (because
        % some of the aa_theoreticalPR values are the same).
        if val1 == sgc                        
            bbcount = bbcount +1;
            % we want to be really sure, so we run the branch and bound
            % algorithm as well (even though that might take some time)
            [status, result2] = system([solverdir,'qapbb < ',filenameinput,' | tail -2']);
            % returns a new permutation as result

            % reformat this new permutation
            [cost2,pos]=textscan(result2,'%d64',1);
            val2=(cons+double(cost2{1}))/N;
            
            if val2 == sgc
                fprintf('fixing (in addition) the following makes the SGC optimal:');
                fixmore(i,:)

                goodcount = goodcount + 1;
                goodfixmore{goodcount}=fixmore(i,:);
            end
        elseif val1 > sgc
            fprintf('this can never happen, the sgc should have been found!')
            return
        end

    end

    save('goodfixings.mat','goodcount','goodfixmore','bbcount','k','ff');
end

return;
