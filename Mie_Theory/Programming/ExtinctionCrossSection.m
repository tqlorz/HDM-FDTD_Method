function C_ext = ExtinctionCrossSection(N, k_o, k_T, k_L, a, epsilon_T, epsilon_m)
    col = length(k_o);
    %% Calculate scatting coefficients
    a_n = zeros(col, 2*N+1);
    for n = -N:N
        a_n(:, n+N+1) = ScattingCoefficients(n, k_o, k_T, k_L, a, epsilon_T, epsilon_m);
    end

    %% Calculate extinction cross-section    
    C_ext = -(2 ./ (k_o * a)) .* sum(real(a_n), 2);
end