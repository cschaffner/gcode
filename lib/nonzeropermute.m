%% permutes non-zero entries of a vector
% takes as input a row vector 
% and returns a matrix whose rows are all possible permutations of the
% non-zero entries of the vector

% by Christian Schaffner, c.schaffner@cwi.nl
% 2 March 2011

function Y = nonzeropermute(X)
  if (size(X,1)>1)
      error('argument needs to be a row vector');
  end
  
  % i is row vector of indices of non-zero positions of X
  % v is row vector of those values
  [ignore,i,v]=find(X);
  
  permmatrix=perms(v);
  
  Y=zeros(size(permmatrix,1),size(X,2));  
  Y(:,i)=permmatrix;
  
end