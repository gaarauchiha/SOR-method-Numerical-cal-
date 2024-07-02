function [x, nit] = sor(A, b, x0, w, d, tol, nmax)
% SOR : solve linear system with SOR iteration
% Usage: [x, nit] = sor(A, b, x0, omega, d, tol, nmax)
% Inputs:
%       A : an n x n-matrix,
%       b : the rhs vector, with length n
%       x0 : the start vector for the iteration
%      tol: error tolerance
%       w: relaxation parameter, (1 < w < 2),
%       d : band width of A.
% Outputs:
%       x : the solution vector
%       nit: number of iterations

n = length(b);
x = x0;
nit = 0;

for k = 1:nmax
    x_old = x;
    for i = 1:n
        sum1 = 0;
        sum2 = 0;
        for j = max(1, i-d):i-1
            sum1 = sum1 + A(i,j) * x(j);
        end
        for j = i+1:min(n, i+d)
            sum2 = sum2 + A(i,j) * x_old(j);
        end
        x(i) = (1-w) * x_old(i) + (w / A(i,i)) * (b(i) - sum1 - sum2); % SOR implementation
    end
    
    % Residual
    r = A * x - b;
    if norm(r, inf) < tol
        break;
    end
    nit = k;
end

if nit == nmax
    warning('Maximum number of iterations reached, Diverges.');
end
end