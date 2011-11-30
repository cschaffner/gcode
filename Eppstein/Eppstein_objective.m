%% performs Eppstein "inverse optimization" for 20 amino properties

% June 29, 2011
% by Christian Schaffner, c.schaffner@uva.nl

function [f g H] = Eppstein_objective(x)

global targetvals

f=-x(21);
for i=1:20
   f=f+ (x(i)-targetvals(i))^2; 
end

if nargout > 1 % gradient required
    g=zeros(21,1);
    for i=1:20
        g(i)=2*(x(i)-targetvals(i));
    end
    g(21)=-1;
    
    if nargout > 2 % Hessian required
        H=zeros(21,21);
        for i=1:20
            H(i,i)=2;
        end
    end
end

end