%% replaces pdist because of stupid license troubles
% does the same as pdist with distance measure of absolute values

% by Christian Schaffner, c.schaffner@uva.nl
% 1 March 2011

function Y = mypdist(X)
%PDIST Pairwise distance between observations.
%   D = PDIST(X) returns a vector D containing the Euclidean distances
%   between each pair of observations in the M-by-N data matrix X. Rows of
%   X correspond to observations, columns correspond to variables. D is a
%   1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%
  n=size(X,1);
  Y=zeros(n,n);
  for i=1:n
      for j=1:n
          Y(i,j)=abs(X(i)-X(j));
      end
  end
  
end