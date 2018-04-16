system("--ticks-per-sec",1000);

int exp7 = 7;
int exp3 = 3;

ring r = 0, (x1, x2, x3, x4, x5, x6, x7), dp;

poly poly1 = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp7) - 1;
poly poly2_1 = (1 + 3*x1*x2 + 5*x2*x3 + 7*x3*x4 + 9*x4*x5 + 11*x5*x6 + 13*x6*x7 + 15*x7*x1)^(exp3);
poly poly2_2 = (1 + 3*x1*x3 + 5*x2*x4 + 7*x3*x5 + 9*x6*x5 + 11*x7*x6 + 13*x6*x1 + 15*x7*x2)^(exp3);
poly poly2_3 = (1 + 3*x1*x4 + 5*x2*x5 + 7*x3*x6 + 9*x6*x7 + 11*x7*x1 + 13*x6*x2 + 15*x7*x3)^(exp3);
poly poly2 = poly2_1 * poly2_2 * poly2_3  - 1;


int t = timer;
list factors = factorize(poly1);
int elapsed = timer - t;
printf("Factor poly1 in Z: %s ms", elapsed);


int t = timer;
list factors = factorize(poly2);
int elapsed = timer - t;
printf("Factor poly2 in Z: %s ms", elapsed);



ring r = 2, (x1, x2, x3, x4, x5, x6, x7), dp;

poly poly1 = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp7) - 1;
poly poly2_1 = (1 + 3*x1*x2 + 5*x2*x3 + 7*x3*x4 + 9*x4*x5 + 11*x5*x6 + 13*x6*x7 + 15*x7*x1)^(exp3);
poly poly2_2 = (1 + 3*x1*x3 + 5*x2*x4 + 7*x3*x5 + 9*x6*x5 + 11*x7*x6 + 13*x6*x1 + 15*x7*x2)^(exp3);
poly poly2_3 = (1 + 3*x1*x4 + 5*x2*x5 + 7*x3*x6 + 9*x6*x7 + 11*x7*x1 + 13*x6*x2 + 15*x7*x3)^(exp3);
poly poly2 = poly2_1 * poly2_2 * poly2_3  - 1;


int t = timer;
list factors = factorize(poly1);
int elapsed = timer - t;
printf("Factor poly1 in Z/2: %s ms", elapsed);


int t = timer;
list factors = factorize(poly2);
int elapsed = timer - t;
printf("Factor poly2 in Z/2: %s ms", elapsed);


ring r = 524287, (x1, x2, x3, x4, x5, x6, x7), dp;

poly poly1 = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^(exp7) - 1;
poly poly2_1 = (1 + 3*x1*x2 + 5*x2*x3 + 7*x3*x4 + 9*x4*x5 + 11*x5*x6 + 13*x6*x7 + 15*x7*x1)^(exp3);
poly poly2_2 = (1 + 3*x1*x3 + 5*x2*x4 + 7*x3*x5 + 9*x6*x5 + 11*x7*x6 + 13*x6*x1 + 15*x7*x2)^(exp3);
poly poly2_3 = (1 + 3*x1*x4 + 5*x2*x5 + 7*x3*x6 + 9*x6*x7 + 11*x7*x1 + 13*x6*x2 + 15*x7*x3)^(exp3);
poly poly2 = poly2_1 * poly2_2 * poly2_3  - 1;


int t = timer;
list factors = factorize(poly1);
int elapsed = timer - t;
printf("Factor poly1 in Z/524287: %s ms", elapsed);


int t = timer;
list factors = factorize(poly2);
int elapsed = timer - t;
printf("Factor poly2 in Z/524287: %s ms", elapsed);



quit;
