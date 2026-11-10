function [limited_solution] = apply_superbee_limiter(solution, dx)
    N = length(solution);
    limited_solution = solution;
    
    for i = 2:N-1
        % Calculate slopes
        left_slope = (solution(i) - solution(i-1)) / dx;
        right_slope = (solution(i+1) - solution(i)) / dx;
        
        % Apply superbee limiter
        limited_slope = superbee_limiter(left_slope, right_slope);
        
        % Update solution with limited slope
        limited_solution(i) = solution(i) + limited_slope * dx/2;
    end
end
