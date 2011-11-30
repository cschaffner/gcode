%% inverts a permutation
% takes as input a row vector with distinct elements (describing a
% permutation)
% and returns the inverse of it

% by Christian Schaffner, c.schaffner@uva.nl
% 24 March 2011

function invp=invertp(p)

  if (size(p,1)>1)
      error('argument needs to be a row vector');
  end

  P=sortrows([p ; 1:size(p,2)]')';
  
  invp=P(2,:);
end
