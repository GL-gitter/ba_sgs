function [A,rhs] = get_SysM2D(p)

[A,rhs] = get_Diff2D(p);
% G_glob = get_Drift2D(p,0);

% A = A+G_glob;

end