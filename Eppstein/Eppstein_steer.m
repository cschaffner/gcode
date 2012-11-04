%% performs Eppstein "inverse optimization" for 20 amino properties

% by Christian Schaffner, c.schaffner@uva.nl
% and Gunnar Klau

% created: May 16, 2011
% cleaned up version: Dec 5, 2011

%% clear workspace 
geneticcode;

%% declare global variables
% in order to use them in Eppstein_objective.m and
% Eppstein_q20_constraint.m
global Bmatrix count popt targetvals

%% set parameters of what we want to do
equif=0;

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

scoretype = 'Eppstein_steer_all1weights';


%% define filenames
solverdir=strcat(fileparts(mfilename('fullpath')),'/../solver/');
Eppsteindir=strcat(fileparts(mfilename('fullpath')),'/');

filenameinput=strcat(solverdir,scoretype, '.input');

filenamegood=strcat(Eppsteindir,scoretype, '_perfect.txt');
filenamepopt=strcat(Eppsteindir,scoretype,'_popt.mat');


%% read in initializing values if they exist

% intialize list of permutations with five permutations that are not the
% identity

popt=zeros(5,20);
popt(1,:)=randperm(20);
popt(2,:)=randperm(20);
popt(3,:)=randperm(20);
popt(4,:)=randperm(20);
popt(5,:)=randperm(20);
while sum(sum((popt==[1:20 ; 1:20 ;1:20 ;1:20 ;1:20]),2)==20)>0
    popt(1,:)=randperm(20);
    popt(2,:)=randperm(20);
    popt(3,:)=randperm(20);
    popt(4,:)=randperm(20);
    popt(5,:)=randperm(20);
end

normTheorPR=[-1.1100   -1.1490   -0.9200   -0.9200   -0.4630    0.0320   -0.5010   -0.4630   -0.3490    0.1090 0.1850    0.5660    0.8320    1.0610    1.8230    2.3560   -1.1870   -0.9580    0.4510    0.6040];


targetvals=normTheorPR;
x=[targetvals';500];
x0=zeros(5,21);
xpolar=targetvals;

goodcount=0;
xgood=zeros(1,21);
pgood=zeros(1,20);
countgood=zeros(1,1);

perfectcount=1;
xperfect=zeros(1,21);
countperfect=zeros(1,1);

if (exist(filenamepopt))
    load(filenamepopt);
    count=size(popt,1)-1;
    popt=popt(1:count,:);
    x0=x0(1:count,:);
    xopt=xopt(1:count,:);
    fprintf('Loaded previous data, starting with %i permutations\n',count)
    x=x0(count,:);
end


%% start the big loop
while 1

count=size(popt,1);
fprintf('calculating with %i permutations now.\n',count)

if (mod(count,10)==0)
    fprintf('reaching %i permutations. saving popt\n',count)
    save(filenamepopt,'popt','xopt','x0','stuckcounter','xvary','perfectcount','xperfect','countperfect');
end

%% do non-convex solver work

% initialization
% "the further we are in the optimization, the more we look into details"
if (x(21) < 50)
    % focus more into the detail by starting at the previous (rounded) end point
    x0(count,:) = x;
else
    % start with the polar requirements to get some broad constraints
    x0(count,:) = [targetvals  500];
end

if (x(21)<5)
  roundconst=10^3;
else 
  roundconst=10^2;
end

% lower bounds
lb =ones(21,1)*(-10);
lb(21) = 0.005;

% upper bounds
ub = ones(21,1)*10;

% equality constraints:

% property values should sum to 0
Aeq = ones(1,21);
Aeq(21)=0;
beq = 0;

% set options
options=optimset('GradConstr','on');
%options.GradObj='on';
options.Jacobian='on';
%options.DerivativeCheck='on';
options.MaxFunEvals=3000;
options.Display='iter';
options.FinDiffType='central';
%options.PlotFcns={@optimplotx,@optimplotfval,@optimplotfirstorderopt};
options.TolFun = 1.000000e-04;
options.TolCon = 1.000000e-04;

[x,fval,exitflag]=fmincon(@Eppstein_objective,x0(count,:),[],[],Aeq,beq,lb,ub,@Eppstein_q20_constraint,options);

% xopt(count) is the value that comes out of the fmincon optimization when
% running starting from x0(count) and 
% using popt(1:count,1:20) as constraints
xopt(count,:)=x;

%% set up problem with new amino properties for solver

% first round x values to three (or two in the beginning) digits after comma:
x=round(x*roundconst)/roundconst;

A = mypdist(x(1:20)') .^ 2;
B = Bmatrix;

CreateQAP(scoretype,[],A,B);

%% run heuristics
% "gqapd.f" GRASP algorithm
[status, result1] = system([solverdir,'gqapd < ',filenameinput,' | tail -2']);
% returns a new permutation as result

% "qapsim.f" from http://www.seas.upenn.edu/qaplib/codes.html
[status, result2] = system([solverdir,'mysimqap < ',filenameinput,' | tail -2']);
% returns a new permutation as result

% reformat this new permutation
[cost1,pos]=textscan(result1,'%d64',1);
pnew1=textscan(result1(pos+1:end),'%2f');
ppnew1=pnew1{1}';

% reformat this new permutation
[cost2,pos]=textscan(result2,'%d64',1);
pnew2=textscan(result2(pos+1:end),'%2f');
ppnew2=pnew2{1}';

% take the best of the two 
if (cost1{1}<=cost2{1})
    ppnew=ppnew1;
    cost=cost1{1};
else
    ppnew=ppnew2;
    cost=cost2{1};
end
    
if (cost<0 || size(ppnew,2)>20 || sum(sort(ppnew) ~= 1:20)>0) 
    fprintf('one of the heuristics provided strange output\n');
    % make a copy of the input file that produced strange input
    system(['cp ',filenameinput,' ',filenameinput,'.',num2str(floor(now*100000))]);
    
    % set the "new" strange permutation to the one previously obtained, so
    % that it does not get added to the list and we only vary the starting values a bit
    ppnew = popt(count,:);
end
    
if (sum(ppnew == 1:20)==20)
    % save these values
    if (perfectcount==0 || (perfectcount > 0 && sum(xperfect(perfectcount,:) == x)) < 20)
        perfectcount=perfectcount+1;        
        xperfect(perfectcount,:)=x;
        countperfect(perfectcount)=count;

        fid = fopen(filenamegood, 'a+');
        fprintf(fid, '\nfound a perfect permutation using %d constraints\n', count);
        fclose(fid);
        dlmwrite(filenamegood, x,'-append','delimiter',' ','precision', '%1.4f');
    end

    fprintf('bingo, not adding this one of course, just varying the values a little and see if we get another bingo.\n');
    ppnew
    stuckcounter=stuckcounter+1;
    xvary(stuckcounter,:)=x;
    if (stuckcounter == 10)
        fprintf('after 10 random tries, we did not find more constraints. We are pretty much done and succeeded.\n');
        xvary
        save(filenamepopt,'popt','xopt','x0','stuckcounter','xvary','perfectcount','xperfect','countperfect');
        break;
    end
    randoffset = (rand(1,20)-0.5)/10;
    x(1:20)=x(1:20)+randoffset;
elseif (sum(ppnew == popt(count,:))<20 && sum(ppnew == popt(count-1,:))<20 && ...
        sum(ppnew == popt(count-2,:))<20 && sum(ppnew == popt(count-3,:))<20 ...
        && sum(ppnew == popt(count-4,:))<20)
    count=count+1;    
    popt(count,:) = ppnew;
    fprintf('added this new permutation\n');
    ppnew
    stuckcounter=0;
    xvary=zeros(20,21);
else
    fprintf('stuck with the following permutation\n');
    ppnew
    % if we are getting back the same permutation
    % try to disturb randomly (by one fourth of a standard deviation, so quite a lot) the starting inputs a bit
    x(1:20)=x(1:20)+(rand(1,20)-0.5)*std(x(1:20))/4;

    stuckcounter=stuckcounter+1;
    xvary(stuckcounter,:)=x;
    if (stuckcounter == 20) 
        fprintf('after 20 random tries, we are pretty much stuck with %i permutations\n',count);
        xvary
        save(filenamepopt,'popt','xopt','x0','stuckcounter','xvary','perfectcount','xperfect','countperfect');
        break;
    end
end


end

