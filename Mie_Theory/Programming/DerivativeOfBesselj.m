function dJ = DerivativeOfBesselj(n, x)
    dJ = 0.5 * (besselj(n - 1, x) - besselj(n + 1, x));
end