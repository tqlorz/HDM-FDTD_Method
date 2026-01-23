function C_ext = ExtinctionCrossSection(k_o, a, a_n)
    C_ext = -(2 ./ (k_o * a)) .* sum(real(a_n), 2);
end