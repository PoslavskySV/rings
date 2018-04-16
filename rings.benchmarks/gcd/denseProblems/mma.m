exp = 7;
a = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp) - 1;
b = (1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7)^(exp) + 1;
g = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7)^(exp) + 3;

ag = Expand[a*g];
bg = Expand[b*g];

r = Timing[PolynomialGCD[ag, bg]];
Print["Done in Z: " <> ToString[ Round[r[[1]]*10^9] ] ];
r = Timing[PolynomialGCD[ag + 1, bg]];
Print["Done in Z (trivial): " <> ToString[ Round[r[[1]]*10^9] ] ];


ag = PolynomialMod[ag, 524287];
bg = PolynomialMod[bg, 524287];

r = Timing[PolynomialGCD[ag, bg, Modulus->524287]];
Print["Done in Zp: " <> ToString[ Round[r[[1]]*10^9] ] ];
r = Timing[PolynomialGCD[ag + 1, bg Modulus->524287]];
Print["Done in Zp (trivial): " <> ToString[ Round[r[[1]]*10^9] ] ];