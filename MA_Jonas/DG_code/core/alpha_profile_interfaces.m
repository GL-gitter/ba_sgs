function alpha = alpha_profile_interfaces(N_chi, alpha_min, nb)
% alpha_profile_interfaces  Symmetric alpha profile on chi-interfaces
% with flat boundary layers (alpha = 1).
%
% alpha = alpha_profile_interfaces(N_chi, alpha_min, nb)
%
% Interfaces: N_chi - 1
%
% - alpha = 1 on the first nb and last nb interfaces
% - smooth (C^1) cosine transition to alpha_min in the center
% - symmetric with respect to the domain center
%
% Inputs:
%   N_chi      number of chi cells (>= 2)
%   alpha_min  minimum alpha in the center (e.g. 0.8)
%   nb         number of boundary interfaces with alpha = 1
%
% Output:
%   alpha      (N_chi-1) x 1 vector

    if nargin < 2 || isempty(alpha_min)
        alpha_min = 0.8;
    end
    if nargin < 3 || isempty(nb)
        nb = 1;
    end

    if N_chi < 2
        error('N_chi must be >= 2.');
    end
    if alpha_min <= 0 || alpha_min > 1
        error('alpha_min must be in (0,1].');
    end

    Nint = N_chi - 1;   % number of interfaces

    % Cap nb so that a transition region remains
    nb = max(0, min(nb, floor((Nint-1)/2)));

    % Initialize
    alpha = ones(Nint,1);

    % Length of transition region on one side
    L = Nint - 2*nb;

    if L <= 0
        % Entire domain is boundary layer
        return;
    end

    % Parametrization of the interior region [0,1]
    s = linspace(0,1,L).';

    % Smooth symmetric dip: 0 at edges of transition, 1 in the middle
    dip = 0.5 * (1 - cos(2*pi*s));   % = sin^2(pi*s)

    % Fill interior
    alpha(nb+1:nb+L) = 1 - (1 - alpha_min) * dip;

    % Safety
    alpha = max(alpha_min, min(1, alpha));
end
