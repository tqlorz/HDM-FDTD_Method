function dH = DerivativeOfHankel1(n, x)
    dH = 0.5 * (besselh(n - 1, 1, x) - besselh(n + 1, 1, x));
end