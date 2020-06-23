

function [w] =  STRidge(X0,y,lam,maxit,tol,normalize)
%     Sequential Threshold Ridge Regression algorithm for finding (hopefully) sparse 
%     approximation to X^{-1}y.  The idea is that this may do better with correlated observables.
%     This assumes y is only one column
if nargin < 6 
   normalize = 2; 
end
d = size(X0,2);
% First normalize data
if normalize ~= 0
    Mreg = zeros(d,1);
    for i = 1:d
        Mreg(i) = 1/(norm(X0(:,i),normalize));
        X(:,i) = Mreg(i)*X0(:,i);
    end
else
    X = X0;
end
   
% Get the standard ridge esitmate
if lam ~= 0
    w = (X'*X + lam*eye(d))\(X'*y);
else
    w =  X\y;
end
num_relevant = d;
biginds = find( abs(w) >= tol);
% Threshold and continue
for j = 1:maxit
    % Figure out which items to cut out
    smallinds = find( abs(w) < tol);
    new_biginds =   setdiff(1:d,smallinds); %[i for i in range(d) if i not in smallinds]
    % If nothing changes then stop
    if num_relevant ==  length(new_biginds);
        break
    else
         num_relevant = length(new_biginds);
    end
    % Also make sure we didn't just lose all the coefficients
    if isempty(new_biginds)
       if j == 1
           return
       else 
           break
       end
    end
    biginds = new_biginds;
    % Otherwise get a new guess
    w(smallinds) = 0;
    if lam ~= 0
         w(biginds) = (X(:, biginds)'*X(:, biginds) + lam*eye(length(biginds)))\(X(:, biginds)'*y);
    else
        w(biginds) = X(:, biginds)\y;
    end
end
% Now that we have the sparsity pattern, use standard least squares to get w
if  ~isempty(biginds)
   w(biginds) = X(:,biginds)\y; 
end
if normalize ~= 0
    w = diag(Mreg)*w;
end


end