%% which assignments need to be fixed in 
% * addition to our 7 such that the SGC is optimal?
% * at all?

% by Christian Schaffner, c.schaffner@uva.nl
% idea by Gunnar Klau

% created: Nov 5, 2012

%% clear workspace 
geneticcode;

%% load results

load goodfixings_all1weights.mat

% potentially fixed blocks
fixed = [1 2 3 10 11 18 19];


for i=1:size(goodfixmore,2)
    if size(intersect(fixed,goodfixmore{i}),2)==size(fixed,2)
        
        i,goodfixmore{i}
    end
end    