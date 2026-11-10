%% Erzeugt das Komplex absorbierende Potential, um Amplituden am Rand zu dämpfen

function W = get_CAP(q, L_xi, delta_xi)
W0 = 1.2;   % 1.2  % 1.4
n = 1.0;    % 1.0  % 2.3
L_xi0 = L_xi/2-delta_xi;
weite = L_xi*0.20;        % 0.2

W = (W0/weite.*(q-L_xi0+weite)).^(2*n).*(q >= L_xi0-weite)...
  + (W0/weite.*(q+L_xi0-weite)).^(2*n).*(q <= -L_xi0+weite);

end
