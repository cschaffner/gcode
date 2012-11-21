function p = randfixperm(n,fixed)
%RANDFIXPERM Random permutation with some positions fixed
%   RANDPERM(n,fixed) is a random permutation of the integers from 1 to n
%   with indices in the array fixed being fixed.
%   For example, RANDPERM(6,[1 3]) might be [1 4 3 5 6 2].
%   
%
% set of indices that are not fixed
notfixed=int8(setdiff(1:n,fixed));

% randomly permute and add a 0 at the beginning
notfixperm=[0 notfixed(randperm(size(notfixed,2)))];

% figure out where the non-fixed positions are
[ignore,locnotfixed]=ismember(1:n,notfixed);

% putting it together again
p = notfixperm(locnotfixed+1)+ int8([1:n] .* ismember(1:n,fixed));

if sum(sort(p) == 1:n)<n
    die('error: randfixperm did not return a permutation');
end

end