package cc.r2.core.poly2;

import cc.r2.core.poly2.DivisionWithRemainder.InverseModMonomial;

import java.util.ArrayList;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.DivisionWithRemainder.quotient;
import static cc.r2.core.poly2.FactorDecomposition.oneFactor;
import static cc.r2.core.poly2.ModularComposition.*;
import static cc.r2.core.poly2.PolynomialArithmetics.polyMultiplyMod;
import static cc.r2.core.poly2.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;


/**
 * Distinct-free factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class DistinctDegreeFactorization {
    private DistinctDegreeFactorization() {}

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} using plain incremental exponents
     * algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationPlain(MutablePolynomialMod poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long factor = poly.lc();
        MutablePolynomialMod base = poly.clone().monic();
        MutablePolynomialMod polyModulus = base.clone();

        if (base.degree <= 1)
            return oneFactor(base, factor);

        if (base.isMonomial())
            return oneFactor(base, factor);

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
        MutablePolynomialMod exponent = MutablePolynomialMod.createMonomial(poly.modulus, 1, 1);
        FactorDecomposition<MutablePolynomialMod> result = new FactorDecomposition<>();
        int i = 0;
        while (!base.isConstant()) {
            ++i;
            exponent = PolynomialArithmetics.polyPowMod(exponent, poly.modulus, polyModulus, invMod, false);
            MutablePolynomialMod tmpExponent = exponent.clone();
            tmpExponent.ensureCapacity(1);
            tmpExponent.data[1] = base.subtractMod(tmpExponent.data[1], 1);
            tmpExponent.fixDegree();
            MutablePolynomialMod gcd = PolynomialGCD(tmpExponent, base);
            if (!gcd.isConstant())
                result.addFactor(gcd.monic(), i);


            base = quotient(base, gcd, false); //can safely destroy reused base
            if (base.degree < 2 * (i + 1)) {// <- early termination
                if (!base.isConstant())
                    result.addFactor(base.monic(), base.degree);
                break;
            }
        }
        return result.setNumericFactor(factor);
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} using
     * plain incremental exponents algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationPrecomputedExponents(MutablePolynomialMod poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long factor = poly.lc();
        MutablePolynomialMod base = poly.clone().monic();
        MutablePolynomialMod polyModulus = base.clone();

        if (base.degree <= 1)
            return oneFactor(base, factor);

        if (base.isMonomial())
            return oneFactor(base, factor);

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
        MutablePolynomialMod exponent = MutablePolynomialMod.createMonomial(poly.modulus, 1, 1);
        FactorDecomposition<MutablePolynomialMod> result = new FactorDecomposition<>();

        ArrayList<MutablePolynomialMod> xPowers = xPowers(polyModulus, invMod);
        int i = 0;
        while (!base.isConstant()) {
            ++i;
            exponent = powModulusMod(exponent, polyModulus, invMod, xPowers);
            MutablePolynomialMod tmpExponent = exponent.clone();
            tmpExponent.ensureCapacity(1);
            tmpExponent.data[1] = poly.subtractMod(tmpExponent.data[1], 1);
            tmpExponent.fixDegree();
            MutablePolynomialMod gcd = PolynomialGCD(tmpExponent, base);
            if (!gcd.isConstant())
                result.addFactor(gcd.monic(), i);

            base = quotient(base, gcd, false); //can safely destroy reused base
            if (base.degree < 2 * (i + 1)) {// <- early termination
                if (!base.isConstant())
                    result.addFactor(base.monic(), base.degree);
                break;
            }
        }
        return result.setNumericFactor(factor);
    }

    /** Baby step / giant step components for d.d.f. in Shoup's algorithm */
    static BabyGiantSteps generateBabyGiantSteps(MutablePolynomialMod poly) {
        int n = poly.degree;
        int B = (int) Math.floor(n / 2.);
        int l = (int) Math.floor(Math.sqrt(B));
        int m = (int) Math.ceil(1.0 * B / l);

        InverseModMonomial invMod = fastDivisionPreConditioning(poly);
        ArrayList<MutablePolynomialMod> xPowers = xPowers(poly, invMod);

        //baby steps
        ArrayList<MutablePolynomialMod> babySteps = new ArrayList<>();
        babySteps.add(MutablePolynomialMod.createMonomial(poly.modulus, 1, 1)); // <- add x
        MutablePolynomialMod xPower = xPowers.get(1); // x^p mod poly
        babySteps.add(xPower); // <- add x^p mod poly
        for (int i = 0; i <= l - 2; ++i)
            babySteps.add(xPower = powModulusMod(xPower, poly, invMod, xPowers));

        // <- xPower = x^(p^l) mod poly

        //giant steps
        ArrayList<MutablePolynomialMod> giantSteps = new ArrayList<>();
        giantSteps.add(MutablePolynomialMod.createMonomial(poly.modulus, 1, 1)); // <- add x
        giantSteps.add(xPower);
        MutablePolynomialMod xPowerBig = xPower;
        int tBrentKung = (int) Math.sqrt(poly.degree);
        ArrayList<MutablePolynomialMod> hPowers = polyPowers(xPowerBig, poly, invMod, tBrentKung);
        for (int i = 0; i < m - 1; ++i)
            giantSteps.add(xPowerBig = compositionBrentKung(xPowerBig, hPowers, poly, invMod, tBrentKung));

        return new BabyGiantSteps(B, l, m, babySteps, giantSteps, invMod);
    }

    /** Shoup's main gcd loop */
    static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationShoup(MutablePolynomialMod poly,
                                                                                      BabyGiantSteps steps) {
        //generate each I_j
        ArrayList<MutablePolynomialMod> iBases = new ArrayList<>();
        for (int j = 0; j <= steps.m; ++j) {
            MutablePolynomialMod iBase = poly.createOne();
            for (int i = 0; i <= steps.l - 1; ++i) {
                MutablePolynomialMod tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i));
                iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, false);
            }
            iBases.add(iBase);
        }

        FactorDecomposition<MutablePolynomialMod> result = new FactorDecomposition<>();

        MutablePolynomialMod current = poly.clone();
        for (int j = 1; j <= steps.m; ++j) {
            MutablePolynomialMod gcd = PolynomialGCD(current, iBases.get(j));
            if (gcd.isConstant())
                continue;
            current = quotient(current, gcd, false);
            for (int i = steps.l - 1; i >= 0; --i) {
                MutablePolynomialMod tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i)));
                if (!tmp.isConstant())
                    result.addFactor(tmp.clone().monic(), steps.l * j - i);

                gcd = quotient(gcd, tmp, false);
            }
        }
        if (!current.isOne())
            result.addFactor(current.monic(), current.degree);

        return result;
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} using Victor Shoup's baby
     * step / giant step algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationShoup(MutablePolynomialMod poly) {
        long factor = poly.lc();
        poly = poly.clone().monic();
        return DistinctDegreeFactorizationShoup(poly, generateBabyGiantSteps(poly)).setNumericFactor(factor);
    }

    /** baby/giant steps for Shoup's d.d.f. algorithm */
    static final class BabyGiantSteps {
        final int B, l, m;
        final ArrayList<MutablePolynomialMod> babySteps;
        final ArrayList<MutablePolynomialMod> giantSteps;
        final InverseModMonomial invMod;

        public BabyGiantSteps(int b, int l, int m, ArrayList<MutablePolynomialMod> babySteps, ArrayList<MutablePolynomialMod> giantSteps, InverseModMonomial invMod) {
            this.B = b;
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly}.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorization(MutablePolynomialMod poly) {
        if (poly.degree < 50)
            return DistinctDegreeFactorizationPlain(poly);
        else if (poly.degree < 200)
            return DistinctDegreeFactorizationPrecomputedExponents(poly);
        else
            return DistinctDegreeFactorizationShoup(poly);
    }

//    static FactorDecomposition<MutablePolynomialMod> fixDistinctDegreeDecomposition(FactorDecomposition<MutablePolynomialMod> f) {
//        for (int i = f.factors.size() - 1; i >= 0; --i) {
//            int exponent = f.exponents.get(i);
//            int degree = f.factors.get(i).degree;
//
//            if (degree % exponent != 0)
//                f.exponents.set(i, degree);
//        }
//        return f;
//    }

    /**
     * Performs square-free factorization followed by distinct-degree factorization modulo {@code modulus}.
     *
     * @param poly the polynomial
     * @return square-free and distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationComplete(MutablePolynomialMod poly) {
        FactorDecomposition<MutablePolynomialMod> squareFree = SquareFreeFactorization(poly);
        long overallFactor = squareFree.factor;

        FactorDecomposition<MutablePolynomialMod> result = new FactorDecomposition<>();
        for (int i = squareFree.factors.size() - 1; i >= 0; --i) {
            FactorDecomposition<MutablePolynomialMod> dd = DistinctDegreeFactorization(squareFree.factors.get(i));
            int nFactors = dd.factors.size();
            for (int j = nFactors - 1; j >= 0; --j)
                result.addFactor(dd.factors.get(j), squareFree.exponents.get(i));
            overallFactor = poly.multiplyMod(overallFactor, dd.factor);
        }

        return result.setNumericFactor(overallFactor);
    }
}
