%% helper function to export QAP problem to solver-readable format

% input: (filename,fixed,A,B)
%   fixed : vector of fixed positions
%   A : 20x20 matrix
%   B : 20x20 matrix
% writes QAP problem to file with name filename in subdirectory 'solver'

% created:    Nov 30, 2011
% by Christian Schaffner, c.schaffner@uva.nl

function CreateQAP(filename,fixed,A,B)

%% define filenames
filenamebase=strcat('solver/',filename);
filenameinput=strcat(filenamebase, '.input');

system(['cd ',pwd]);
%filenameout=strcat('solver/',scoretype, '_'.output');


%% export matrix to solver

nrfixed=size(fixed,2);

% initialize empty list of indices to permute
indper=[];
for i=1:20
    if (sum(ismember(i, fixed))==0) 
        indper=[indper,i];
    end
end

NA=1;
NB=1;

% convert to integers
while (norm(A*NA-round(A*NA)) > 0.0001)  NA=NA*10; end
fprintf('A matrix was multiplied by: %d \n',NA);
if (NA>10^6) % then A should be rounded to start with!
    error('A should be rounded first, solver does not like HUGE integers');
end

while (norm(B*NB-round(B*NB)) > 0.0001)  NB=NB*10; end
fprintf('B matrix was multiplied by: %d \n',NB);
if (NB>10^6) % then B should be rounded to start with!
    error('B should be rounded first, solver does not like HUGE integers');
end

% restricted matrices to export
Aexp=A(indper,indper)*NA;
Bexp=B(indper,indper)*NB;

if (min(min(Aexp))<0)
    minn= -min(min(Aexp));
    % make all entries positive
    Aexp = Aexp + minn;
    fprintf('%d has been added to all entries of the A matrix \n',minn);
    fprintf('this change increases the objective value by %d \n',minn*sum(sum(B)));    
end

if (min(min(Bexp))<0)
    error('the B matrix should not contain negative values...');
end

% constant term
fprintf('constant term: %6.0f \n',sum(sum(A(fixed,fixed) .* B(fixed,fixed)*NA*NB)))

% linear term (times two, because symmetric)
Cexp=2*A(sort(fixed),indper)'*B(sort(fixed),indper)*NA*NB;

n=size(Aexp,1);
dlmwrite(filenameinput, n, 'delimiter', ' ', 'precision', '%6.0f');
dlmwrite(filenameinput, Aexp, '-append', 'delimiter', ' ', 'precision', '%6.0f', 'roffset', 1);
dlmwrite(filenameinput, Bexp, '-append', 'delimiter', ' ', 'precision', '%6.0f', 'roffset', 1);
if (nrfixed>0)
    dlmwrite(filenameinput, Cexp, '-append', 'delimiter', ' ', 'precision', '%6.0f', 'roffset', 1);    
end

end
