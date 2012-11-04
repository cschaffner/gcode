%% takes a matrix and normalizes the rows of it
% by making the mean zero and the standard deviation 1

% by Christian Schaffner, c.schaffner@uva.nl
% June 9, 2011

function Y = normalize(X)
  n=size(X,1);
  if (size(X,2) ~= 20)
      fprintf('for amino acids, we are expecting a matrix with 20 columns');
  end
  Y=zeros(n,size(X,2));
  for i=1:n
      Y(i,:) = (X(i,:) - mean(X(i,:)) ) / std(X(i,:));
  end  
end