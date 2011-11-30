% takes a matrix of values where columns stand for data properties
% and a row vector of according weights
% returns a matrix with weighted Euclidean distances between pairs of
% values as in Higgs 2009

% d(a,b) = (sum_k w_k (p_ka - p_kb)^2 )^(1/2)

% by Christian Schaffner, c.schaffner@cwi.nl
% 1 March 2011

function Y = mypdistweights(X,w)

  n=size(X,1);
  m=size(X,2);
  
  if (m ~= size(w,2))
      error('number of weights has to match number of columns of the matrix');
  end
  
  Y=zeros(n,n);
  for i=1:n
      for j=1:n
          for k=1:m
             Y(i,j) = Y(i,j) + (w(k) * (X(i,k)-X(j,k))^2 ) ;
          end
          Y(i,j) =  sqrt(Y(i,j));
      end
  end
  
end