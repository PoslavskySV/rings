package cc.r2.core.poly2;

import static cc.r2.core.poly2.DistinctDegreeFactorization.DistinctDegreeFactorization;
import static cc.r2.core.poly2.EqualDegreeFactorization.CantorZassenhaus;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;

/**
 * Factorization of univariate polynomials with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Factorization {
    private Factorization() {}


    /** x^n * poly */
    private static final class FactorMonomial<T extends MutablePolynomialAbstract<T>> {
        final T theRest, monomial;

        FactorMonomial(T theRest, T monomial) {
            this.theRest = theRest;
            this.monomial = monomial;
        }
    }

    /** factor out common monomial term (x^n) */
    private static <T extends MutablePolynomialAbstract<T>> FactorMonomial<T> factorOutMonomial(T poly) {
        int i = 0;
        while (poly.data[i] == 0) ++i;
        assert i < poly.data.length;

        if (i == 0)
            return new FactorMonomial<>(poly, poly.createOne());
        return new FactorMonomial<>(poly.clone().shiftLeft(i), poly.createMonomial(1, i));
    }


    /** early check for trivial cases */
    private static <T extends MutablePolynomialAbstract<T>> FactorDecomposition<T> earlyFactorizationChecks(T poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return FactorDecomposition.oneFactor(poly, 1);

        return null;
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static FactorDecomposition<MutablePolynomialMod> factor(MutablePolynomialMod poly) {
        FactorDecomposition<MutablePolynomialMod> result = earlyFactorizationChecks(poly);
        if (result != null)
            return result;

        FactorMonomial<MutablePolynomialMod> base = factorOutMonomial(poly);

        result = new FactorDecomposition<>();
        result.addFactor(base.monomial, 1);

        //do square-free factorization
        FactorDecomposition<MutablePolynomialMod> sqf = SquareFreeFactorization(base.theRest);
        for (int i = 0; i < sqf.factors.size(); ++i) {
            //for each square-free factor
            MutablePolynomialMod sqfFactor = sqf.factors.get(i);
            int sqfExponent = sqf.exponents.get(i);

            //do distinct-degree factorization
            FactorDecomposition<MutablePolynomialMod> ddf = DistinctDegreeFactorization(sqfFactor);
            for (int j = 0; j < ddf.factors.size(); ++j) {
                //for each distinct-degree factor
                MutablePolynomialMod ddfFactor = ddf.factors.get(j);
                int ddfExponent = ddf.exponents.get(j);

                //do equal-degree factorization
                FactorDecomposition<MutablePolynomialMod> edf = CantorZassenhaus(ddfFactor, ddfExponent);
                for (MutablePolynomialMod irreducibleFactor : edf.factors)
                    //put final irreducible factor into the result
                    result.addFactor(irreducibleFactor.monic(), sqfExponent);
            }
        }

        return result.setNumericFactor(poly.lc());
    }
}
