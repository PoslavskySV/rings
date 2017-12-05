package cc.redberry.rings.poly.univar;

/**
 * @author Stanislav Poslavsky
 * @since 2.1
 */
public final class DiophantineEquations {
    private DiophantineEquations() {}

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    public static <Poly extends IUnivariatePolynomial<Poly>>
    Poly[] monicExtendedEuclid(Poly a, Poly b) {
        Poly[] xgcd = UnivariateGCD.PolynomialExtendedGCD(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        //normalize: x * a + y * b = 1
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[0].monic();

        return xgcd;
    }

    /**
     * Solves a1 * x1 + a2 * x2 + ... = rhs for given univariate and rhs and unknown x_i
     */
    public static final class DiophantineSolver<Poly extends IUnivariatePolynomial<Poly>> {
        /** the given factors */
        final Poly[] factors;
        final Poly[] solution;
        final Poly gcd;

        @SuppressWarnings("unchecked")
        public DiophantineSolver(Poly[] factors) {
            this.factors = factors;
            this.solution = factors[0].createArray(factors.length);

            Poly prev = factors[0];
            solution[0] = factors[0].createOne();

            for (int i = 1; i < factors.length; i++) {
                Poly[] xgcd = monicExtendedEuclid(prev, factors[i]);
                for (int j = 0; j < i; j++)
                    solution[j].multiply(xgcd[1]);
                solution[i] = xgcd[2];
                prev = xgcd[0];
            }
            gcd = prev;
        }

        public Poly[] solve(Poly rhs) {
            rhs = UnivariateDivision.divideOrNull(rhs, gcd, true);
            if (rhs == null)
                throw new IllegalArgumentException("Not solvable.");
            Poly[] solution = rhs.createArray(this.solution.length);
            for (int i = 0; i < solution.length; i++)
                solution[i] = this.solution[i].clone().multiply(rhs);
            return solution;
        }
    }
}
