%% boolean function on inputs 1<=i,j<=4
% returns 1 if the change from i to j is a transition (U-C,C-U or A-G,G-A)
% and 0 if it's a transversion

%  where U=1, C=2, A=3, G=4

% by Christian Schaffner, c.schaffner@cwi.nl
% 28 February 2011

function transit=istransit(i,j)

  transit=0;
  
  if (i>j)
      k=i; i=j; j=k; % swap i and j
  elseif(i==j)
      error('inputs to istransit function have to differ');
  end
      
  if (abs(i-j)==1 && (i==1 || i==3) )
      transit=1;
  end
end
