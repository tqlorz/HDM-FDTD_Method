function a_n = ScattingCoefficients_Nonlocal(N, k_o, k_T, k_L, a, epsilon_T, epsilon_m)
    raw = length(k_o);
    a_n = zeros(raw, 2*N+1);

    for n = -N:N
        % Calculate Bessel functions and their derivatives
        J_n_koa = besselj(n, k_o * a);
        J_n_kTa = besselj(n, k_T * a);
        J_n_kLa = besselj(n, k_L * a);
        H_n_koa = besselh(n, 1, k_o * a);
        dJ_n_koa = DerivativeOfBesselj(n, k_o * a);
        dJ_n_kTa = DerivativeOfBesselj(n, k_T * a);
        dJ_n_kLa = DerivativeOfBesselj(n, k_L * a);
        dH_n_koa = DerivativeOfHankel1(n, k_o * a);
    
        % Calculate auxiliary valriables c_n
        c_n = n^2./(k_L*a) .* J_n_kLa./dJ_n_kLa .* J_n_kTa ...
                            .* (sqrt(epsilon_T)./(k_o*a) - sqrt(epsilon_m)./(k_T*a));
    
        % Calculate the scattering coefficient a_n
        a_n(:, n+N+1) = -(sqrt(epsilon_m) .* J_n_koa .* dJ_n_kTa - sqrt(epsilon_T) .* dJ_n_koa .* J_n_kTa ...
                + c_n .* J_n_koa) ./ ...
              (sqrt(epsilon_m) .* H_n_koa .* dJ_n_kTa - sqrt(epsilon_T) .* dH_n_koa .* J_n_kTa ...
                + c_n .* H_n_koa);
    end
end