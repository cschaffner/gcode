%% check which fixed amino-acid assignments make the standard genetic code optimal
% when using the subset approach

% file expects as input the output of permutecode_allsubsets.m which checks
% all possible 120^4 permutations (when using subsets) and stores the ones
% that are better than the SGC

% created: March 7, 2011
% cleaned up: Nov 3, 2012

% by Christian Schaffner, c.schaffner@uva.nl

%% initialize stuff
geneticcode;

% potentially fixed blocks
%fixed = [1 2 3 10 11 14 18 19];
fixed = [1 2 3 10 11 18 19];


% load results of all subset computations
%scoretype = 'theoretic Polar all subsets';
%load 'th_polar_all1weights_all_good_subset_permutations.mat'
load 'th_polar_FH_weights_all_good_subset_permutations.mat'

% scoretype = 'Benner all subsets';
% load 'output/Benner all subsets0blcksfix1samples.mat'

% reads in A,B,Bnorm,count,minper,minvals,means

% get rids of tons of zeros
minper=minper(1:count,:);
minvals=minvals(1:count,:);


% standard weights:
wtransit1=1;
wtransver1=1;
wtransit2=1;
wtransver2=1;
wtransit3=1;
wtransver3=1;

%% do the work

disp('2 blocks fixed');
% create all possibilities of 2 fixed blocks
fixblocks=mycombnk(fixed,2);
for i=1:size(fixblocks,1)
    fb=fixblocks(i,:);
    % filter out the codes that are better than the SGC for fixed blocks fb
    bettercodes=find(minper(:,fb(1))==fb(1) & minper(:,fb(2))==fb(2));
    nb=size(bettercodes,1);
    if (nb<5 && min(minvals(bettercodes))==sgc(1) )
        fb   
        bettercodes
        minvals(bettercodes)
    end
end

disp('3 blocks fixed');
% create all possibilities of 3 fixed blocks
fixblocks=mycombnk(fixed,3);
for i=1:size(fixblocks,1)
    fb=fixblocks(i,:);
    % filter out the codes that are better than the SGC for fixed blocks fb
    bettercodes=find(minper(:,fb(1))==fb(1) & minper(:,fb(2))==fb(2) & minper(:,fb(3))==fb(3));
    nb=size(bettercodes,1);
    if (nb<5 && min(minvals(bettercodes))==sgc(1))
        fb   
        bettercodes
        minvals(bettercodes)
    end
end

disp('4 blocks fixed');
fixblocks=mycombnk(fixed,4);
for i=1:size(fixblocks,1)
    fb=fixblocks(i,:);
    % filter out the codes that are better than the SGC for fixed blocks fb
    bettercodes=find(minper(:,fb(1))==fb(1) & minper(:,fb(2))==fb(2)& minper(:,fb(3))==fb(3) & minper(:,fb(4))==fb(4));
    nb=size(bettercodes,1);
    if (nb<5 && min(minvals(bettercodes))==sgc(1))
        fb   
        bettercodes
        minvals(bettercodes)
    end
end
% 
% fixblocks=mycombnk(fixed,5);
% for i=1:size(fixblocks,1)
%     fb=fixblocks(i,:);
%     nb=size(find(minper(:,fb(1))==fb(1) & minper(:,fb(2))==fb(2) & minper(:,fb(3))==fb(3) & minper(:,fb(4))==fb(4) & minper(:,fb(5))==fb(5)),1);
%     if (nb==1)
%         fb        
%     end
% end