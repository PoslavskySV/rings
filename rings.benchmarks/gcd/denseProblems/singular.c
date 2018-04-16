system("--ticks-per-sec",1000);
ring r = 0, (x1, x2, x3, x4, x5, x6, x7), dp;

int exp = 7;
poly a = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp) - 1;
poly b = (1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7)^(exp) + 1;
poly g = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7)^(exp) + 3;

poly ag = a*g;
poly bg = b*g;

int t = timer;
poly gcdres = gcd(ag, bg);
int elapsed = timer - t;
printf("Done in Z: %s", elapsed);


int t = timer;
poly gcdres = gcd(ag + 1, bg);
int elapsed = timer - t;
printf("Done in Z (trivial): %s", elapsed);


ring r = 524287, (x1, x2, x3, x4, x5, x6, x7), dp;

poly a = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp) - 1;
poly b = (1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7)^(exp) + 1;
poly g = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7)^(exp) + 3;

poly ag = a*g;
poly bg = b*g;


int t = timer;
poly gcdres = gcd(ag, bg);
int elapsed = timer - t;
printf("Done in Zp: %s", elapsed);


int t = timer;
poly gcdres = gcd(ag + 1, bg);
int elapsed = timer - t;
printf("Done in Zp (trivial): %s", elapsed);

quit;
