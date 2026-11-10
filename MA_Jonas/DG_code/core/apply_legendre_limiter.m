function [limited_solution] = apply_legendre_limiter(solution, k, alpha)
% function [limited_solution] = apply_legendre_limiter(solution, k, alpha)
% Purpose: Apply the Legendre limiter  to a solution vector.
% Variables: solution vector, nodes/element, numerical viscosity factor
% Notes: high-order
% Paper: https://doi.org/10.1016/j.jcp.2007.05.011
limited_solution = solution;
r = linspace(-1,1,k);
N = length(solution);
L = zeros(k, k);
c = zeros(N, 1);
for i = 1:k
    L(:,i) = legendreP(i-1, r);
end

for i = 1:k:N
    rhs = solution(i:i+k-1);
    c_element = L \ rhs;
    c(i:i+k-1) = flip(c_element);
end

for s = 2:N/k-1
    j = (s-1)*k+1;

    for i = j:j+k-2 

        c_neu = minmod_limiter(c(i), alpha*(c(i+k+1)-c(i+1)), alpha*(c(i+1)-c(i-k+1)));
        
        if c(i) == c_neu
            break;
        
        else 
            c(i) = c_neu;
        end
 
    end

    limited_solution(j:j+k-1) = L*flip(c(j:j+k-1));

end

end

