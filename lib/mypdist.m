%% replaces pdist because of stupid license troubles
% does the same as pdist with distance measure of absolute values

% by Christian Schaffner, c.schaffner@uva.nl
% 1 March 2011

function Y = mypdist(X)
%PDIST Pairwise distance between observations.

  n=size(X,1);
  Y=zeros(n,n);
  for i=1:n
      for j=1:n
          Y(i,j)=abs(X(i)-X(j));
      end
  end
  
end