function [limited_solution] = apply_minmod_limiter(solution, dx)
    N = length(solution);
    limited_solution = solution;
    
    for i = 2:N-1
        % Calculate slopes
        left_slope = (solution(i) - solution(i-1)) / dx;
        right_slope = (solution(i+1) - solution(i)) / dx;
        center_slope = (solution(i+1) - solution(i-1)) / (2*dx);
        
        % Apply minmod limiter
        limited_slope = minmod_limiter(left_slope, right_slope, center_slope);
        
        % Update solution with limited slope
        limited_solution(i) = solution(i) + limited_slope * dx/2;
    end
end