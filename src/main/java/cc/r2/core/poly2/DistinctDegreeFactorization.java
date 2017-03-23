package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly2.DivisionWithRemainder.InverseModMonomial;

import java.util.ArrayList;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.DivisionWithRemainder.quotient;
import static cc.r2.core.poly2.ModularComposition.*;
import static cc.r2.core.poly2.PolynomialArithmetics.polyMultiplyMod;
import static cc.r2.core.poly2.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;
import static cc.r2.core.poly2.lFactorDecomposition.oneFactor;


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
    public static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationPlain(MutablePolynomialMod poly) {
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
        lFactorDecomposition<MutablePolynomialMod> result = new lFactorDecomposition<>();
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
    public static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationPrecomputedExponents(MutablePolynomialMod poly) {
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
        lFactorDecomposition<MutablePolynomialMod> result = new lFactorDecomposition<>();

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

    /** Shoup's parameter */
    private static final double SHOUP_BETA = 0.5;

//    /** Baby step / giant step components for d.d.f. in Shoup's algorithm */
//    static BabyGiantSteps generateBabyGiantSteps(MutablePolynomialMod poly) {
//        int n = poly.degree;
//        int l = (int) Math.ceil(Math.pow(1.0 * n, SHOUP_BETA));
//        int m = (int) Math.ceil(1.0 * n / 2 / l);
//
//        InverseModMonomial invMod = fastDivisionPreConditioning(poly);
//        ArrayList<MutablePolynomialMod> xPowers = xPowers(poly, invMod);
//
//        //baby steps
//        ArrayList<MutablePolynomialMod> babySteps = new ArrayList<>();
//        babySteps.add(MutablePolynomialMod.createMonomial(poly.modulus, 1, 1)); // <- add x
//        MutablePolynomialMod xPower = xPowers.get(1); // x^p mod poly
//        babySteps.add(xPower); // <- add x^p mod poly
//        for (int i = 0; i <= l - 2; ++i)
//            babySteps.add(xPower = powModulusMod(xPower, poly, invMod, xPowers));
//
//        // <- xPower = x^(p^l) mod poly
//
//        //giant steps
//        ArrayList<MutablePolynomialMod> giantSteps = new ArrayList<>();
//        giantSteps.add(MutablePolynomialMod.createMonomial(poly.modulus, 1, 1)); // <- add x
//        giantSteps.add(xPower);
//        MutablePolynomialMod xPowerBig = xPower;
//        int tBrentKung = (int) Math.sqrt(poly.degree);
//        ArrayList<MutablePolynomialMod> hPowers = polyPowers(xPowerBig, poly, invMod, tBrentKung);
//        for (int i = 0; i < m - 1; ++i)
//            giantSteps.add(xPowerBig = compositionBrentKung(xPowerBig, hPowers, poly, invMod, tBrentKung));
//
//        return new BabyGiantSteps(l, m, babySteps, giantSteps, invMod);
//    }

    /** Baby step / giant step components for d.d.f. in Shoup's algorithm */
    static <T extends IMutablePolynomialZp<T>> BabyGiantSteps<T> generateBabyGiantSteps(T poly) {
        int n = poly.degree();
        int l = (int) Math.ceil(Math.pow(1.0 * n, SHOUP_BETA));
        int m = (int) Math.ceil(1.0 * n / 2 / l);

        InverseModMonomial<T> invMod = fastDivisionPreConditioning(poly);
        ArrayList<T> xPowers = xPowers(poly, invMod);

        //baby steps
        ArrayList<T> babySteps = new ArrayList<>();
        babySteps.add(poly.createMonomial(1)); // <- add x
        T xPower = xPowers.get(1); // x^p mod poly
        babySteps.add(xPower); // <- add x^p mod poly
        for (int i = 0; i <= l - 2; ++i)
            babySteps.add(xPower = powModulusMod(xPower, poly, invMod, xPowers));

        // <- xPower = x^(p^l) mod poly

        //giant steps
        ArrayList<T> giantSteps = new ArrayList<>();
        giantSteps.add(poly.createMonomial(1)); // <- add x
        giantSteps.add(xPower);
        T xPowerBig = xPower;
        int tBrentKung = (int) Math.sqrt(poly.degree());
        ArrayList<T> hPowers = polyPowers(xPowerBig, poly, invMod, tBrentKung);
        for (int i = 0; i < m - 1; ++i)
            giantSteps.add(xPowerBig = compositionBrentKung(xPowerBig, hPowers, poly, invMod, tBrentKung));

        return new BabyGiantSteps<>(l, m, babySteps, giantSteps, invMod);
    }

//    /** Shoup's main gcd loop */
//    static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationShoup(MutablePolynomialMod poly,
//                                                                                      BabyGiantSteps steps) {
//        //generate each I_j
//        ArrayList<MutablePolynomialMod> iBases = new ArrayList<>();
//        for (int j = 0; j <= steps.m; ++j) {
//            MutablePolynomialMod iBase = poly.createOne();
//            for (int i = 0; i <= steps.l - 1; ++i) {
//                MutablePolynomialMod tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i));
//                iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, false);
//            }
//            iBases.add(iBase);
//        }
//
//        lFactorDecomposition<MutablePolynomialMod> result = new lFactorDecomposition<>();
//
//        MutablePolynomialMod current = poly.clone();
//        for (int j = 1; j <= steps.m; ++j) {
//            MutablePolynomialMod gcd = PolynomialGCD(current, iBases.get(j));
//            if (gcd.isConstant())
//                continue;
//            current = quotient(current, gcd, false);
//            for (int i = steps.l - 1; i >= 0; --i) {
//                MutablePolynomialMod tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i)));
//                if (!tmp.isConstant())
//                    result.addFactor(tmp.clone().monic(), steps.l * j - i);
//
//                gcd = quotient(gcd, tmp, false);
//            }
//        }
//        if (!current.isOne())
//            result.addFactor(current.monic(), current.degree);
//
//        return result;
//    }

    /** Shoup's main gcd loop */
    static <T extends IMutablePolynomialZp<T>> void DistinctDegreeFactorizationShoup(T poly,
                                                                                     BabyGiantSteps<T> steps,
                                                                                     FactorDecomposition<T> result) {
        //generate each I_j
        ArrayList<T> iBases = new ArrayList<>();
        for (int j = 0; j <= steps.m; ++j) {
            T iBase = poly.createOne();
            for (int i = 0; i <= steps.l - 1; ++i) {
                T tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i));
                iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, false);
            }
            iBases.add(iBase);
        }

        T current = poly.clone();
        for (int j = 1; j <= steps.m; ++j) {
            T gcd = PolynomialGCD(current, iBases.get(j));
            if (gcd.isConstant())
                continue;
            current = quotient(current, gcd, false);
            for (int i = steps.l - 1; i >= 0; --i) {
                T tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i)));
                if (!tmp.isConstant())
                    result.addFactor(tmp.clone().monic(), steps.l * j - i);

                gcd = quotient(gcd, tmp, false);
            }
        }
        if (!current.isOne())
            result.addFactor(current.monic(), current.degree());
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
    public static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationShoup(MutablePolynomialMod poly) {
        long factor = poly.lc();
        poly = poly.clone().monic();
        lFactorDecomposition<MutablePolynomialMod> result = new lFactorDecomposition<>();
        DistinctDegreeFactorizationShoup(poly, generateBabyGiantSteps(poly), result);
        return result.setNumericFactor(factor);
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
    public static bFactorDecomposition<bMutablePolynomialMod> DistinctDegreeFactorizationShoup(bMutablePolynomialMod poly) {
        BigInteger factor = poly.lc();
        poly = poly.clone().monic();
        bFactorDecomposition<bMutablePolynomialMod> result = new bFactorDecomposition<>();
        DistinctDegreeFactorizationShoup(poly, generateBabyGiantSteps(poly), result);
        return result.setNumericFactor(factor);
    }

    /** baby/giant steps for Shoup's d.d.f. algorithm */
    static final class BabyGiantSteps<T extends IMutablePolynomialZp<T>> {
        final int l, m;
        final ArrayList<T> babySteps;
        final ArrayList<T> giantSteps;
        final InverseModMonomial<T> invMod;

        public BabyGiantSteps(int l, int m, ArrayList<T> babySteps, ArrayList<T> giantSteps, InverseModMonomial<T> invMod) {
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    /** when to switch to Shoup's algorithm */
    private static final int DEGREE_SWITCH_TO_SHOUP = 256;

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly}.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorization(MutablePolynomialMod poly) {
        if (poly.degree < DEGREE_SWITCH_TO_SHOUP)
            return DistinctDegreeFactorizationPrecomputedExponents(poly);
        else
            return DistinctDegreeFactorizationShoup(poly);
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly}.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    public static bFactorDecomposition<bMutablePolynomialMod> DistinctDegreeFactorization(bMutablePolynomialMod poly) {
        return DistinctDegreeFactorizationShoup(poly);
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly}.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may be incomplete.
     *
     * @param poly the polynomial
     * @return distinct-degree decomposition of {@code poly}
     */
    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomialZp<T>> FactorDecomposition<T> DistinctDegreeFactorization(T poly) {
        if (poly instanceof MutablePolynomialMod)
            return (FactorDecomposition<T>) DistinctDegreeFactorization((MutablePolynomialMod) poly);
        else
            return (FactorDecomposition<T>) DistinctDegreeFactorization((bMutablePolynomialMod) poly);
    }

//    static lFactorDecomposition<MutablePolynomialMod> fixDistinctDegreeDecomposition(lFactorDecomposition<MutablePolynomialMod> f) {
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
    public static lFactorDecomposition<MutablePolynomialMod> DistinctDegreeFactorizationComplete(MutablePolynomialMod poly) {
        lFactorDecomposition<MutablePolynomialMod> squareFree = SquareFreeFactorization(poly);
        long overallFactor = squareFree.factor;

        lFactorDecomposition<MutablePolynomialMod> result = new lFactorDecomposition<>();
        for (int i = squareFree.size() - 1; i >= 0; --i) {
            lFactorDecomposition<MutablePolynomialMod> dd = DistinctDegreeFactorization(squareFree.get(i));
            int nFactors = dd.size();
            for (int j = nFactors - 1; j >= 0; --j)
                result.addFactor(dd.get(j), squareFree.getExponent(i));
            overallFactor = poly.multiplyMod(overallFactor, dd.factor);
        }

        return result.setNumericFactor(overallFactor);
    }
}
