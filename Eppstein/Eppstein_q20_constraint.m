%% performs Eppstein "inverse optimization" for 20 amino properties

% created: May 16, 2011
% by Christian Schaffner, c.schaffner@uva.nl
% and Gunnar Klau

% cleaned up on: June 8, 2011

function [c, ceq,gradc,gradceq] = Eppstein_q20_constraint(x)

global Bmatrix popt count

Y=zeros(20,20);
for i=1:20
  for j=1:20
      Y(i,j)=(x(i)-x(j))^2;
  end
end

% compute constraint values
c = zeros(count,1);
for l=1:count
  c(l) = - ( sum(sum( Y .* (Bmatrix(popt(l,:),popt(l,:))-Bmatrix) )) - x(21)  ) ;
end

ceq = [];

% compute gradient constraint values
if nargout > 2
    gradc=zeros(21,count);
    for l=1:count
        for i=1:20 % derivation with respect to x(i)
            Y=zeros(20,20);
            for j=1:20
                Y(i,j)=2*(x(i)-x(j));
                Y(j,i)=2*(x(i)-x(j));
            end
            % derivation of c(l) with respect to x(i)
            gradc(i,l)=- (sum(sum( Y .* (Bmatrix(popt(l,:),popt(l,:))-Bmatrix) )) );
        end
        gradc(21,l) = 1;
    end
    gradceq = [];
end