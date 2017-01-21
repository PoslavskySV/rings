package cc.r2.core.polynomial;

import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;

import static cc.r2.core.polynomial.DivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.DivideAndRemainder.quotient;
import static cc.r2.core.polynomial.FactorDecomposition.oneFactor;
import static cc.r2.core.polynomial.LongArithmetics.*;
import static cc.r2.core.polynomial.ModularComposition.*;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyMultiplyMod;
import static cc.r2.core.polynomial.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.polynomial.SquareFreeFactorization.SquareFreeFactorization;

/**
 * Distinct-free factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 */
public final class DistinctDegreeFactorization {
    private DistinctDegreeFactorization() {}

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} modulo {@code modulus} using
     * plain incremental exponents algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may not be complete.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    public static FactorDecomposition DistinctDegreeFactorizationPlain(
            MutablePolynomial poly, long modulus) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long factor = mod(poly.lc(), modulus);
        MutablePolynomial base = poly.clone().monic(modulus);
        MutablePolynomial polyModulus = base.clone();

        if (base.degree <= 1)
            return oneFactor(base, factor);

        if (base.isMonomial())
            return oneFactor(base, factor);

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        MutablePolynomial exponent = MutablePolynomial.create(0, 1);
        ArrayList<MutablePolynomial> factors = new ArrayList<>();
        TIntArrayList degrees = new TIntArrayList();
        int i = 0;
        while (!base.isConstant()) {
            ++i;
            exponent = PolynomialArithmetics.polyPowMod(exponent, modulus, polyModulus, invMod, modulus, false);
            MutablePolynomial tmpExponent = exponent.clone();
            tmpExponent.ensureCapacity(1);
            tmpExponent.data[1] = subtractMod(tmpExponent.data[1], 1, modulus);
            tmpExponent.fixDegree();
            MutablePolynomial gcd = PolynomialGCD(tmpExponent, base, modulus);
            if (!gcd.isConstant()) {
                factors.add(gcd.monic(modulus));
                degrees.add(i);
            }
            base = quotient(base, gcd, modulus, false); //can safely destroy reused base
            if (base.degree < 2 * (i + 1)) {// <- early termination
                if (!base.isConstant()) {
                    factors.add(base.monic(modulus));
                    degrees.add(base.degree);
                }
                break;
            }
        }
        return new FactorDecomposition(factors, degrees, factor);
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} modulo {@code modulus} using
     * plain incremental exponents algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may not be complete.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    public static FactorDecomposition DistinctDegreeFactorizationPrecomputedExponents(
            MutablePolynomial poly, long modulus) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long factor = mod(poly.lc(), modulus);
        MutablePolynomial base = poly.clone().monic(modulus);
        MutablePolynomial polyModulus = base.clone();

        if (base.degree <= 1)
            return oneFactor(base, factor);

        if (base.isMonomial())
            return oneFactor(base, factor);

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        MutablePolynomial exponent = MutablePolynomial.create(0, 1);
        ArrayList<MutablePolynomial> factors = new ArrayList<>();
        TIntArrayList degrees = new TIntArrayList();

        ArrayList<MutablePolynomial> xPowers = xPowers(polyModulus, invMod, modulus);
        int i = 0;
        while (!base.isConstant()) {
            ++i;
            exponent = powModulusMod(exponent, polyModulus, invMod, modulus, xPowers);
            MutablePolynomial tmpExponent = exponent.clone();
            tmpExponent.ensureCapacity(1);
            tmpExponent.data[1] = subtractMod(tmpExponent.data[1], 1, modulus);
            tmpExponent.fixDegree();
            MutablePolynomial gcd = PolynomialGCD(tmpExponent, base, modulus);
            if (!gcd.isConstant()) {
                factors.add(gcd.monic(modulus));
                degrees.add(i);
            }
            base = quotient(base, gcd, modulus, false); //can safely destroy reused base
            if (base.degree < 2 * (i + 1)) {// <- early termination
                if (!base.isConstant()) {
                    factors.add(base.monic(modulus));
                    degrees.add(base.degree);
                }
                break;
            }
        }
        return new FactorDecomposition(factors, degrees, factor);
    }

    /** Baby step / giant step components for d.d.f. in Shoup's algorithm */
    static BabyGiantSteps generateBabyGiantSteps(MutablePolynomial poly, long modulus) {
        int n = poly.degree;
        int B = (int) Math.floor(n / 2.);
        int l = (int) Math.floor(Math.sqrt(B));
        int m = (int) Math.ceil(1.0 * B / l);

        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        ArrayList<MutablePolynomial> xPowers = xPowers(poly, invMod, modulus);

        //baby steps
        ArrayList<MutablePolynomial> babySteps = new ArrayList<>();
        babySteps.add(MutablePolynomial.createMonomial(1, 1)); // <- add x
        MutablePolynomial xPower = xPowers.get(1); // x^p mod poly
        babySteps.add(xPower); // <- add x^p mod poly
        for (int i = 0; i <= l - 2; ++i)
            babySteps.add(xPower = powModulusMod(xPower, poly, invMod, modulus, xPowers));

        // <- xPower = x^(p^l) mod poly

        //giant steps
        ArrayList<MutablePolynomial> giantSteps = new ArrayList<>();
        giantSteps.add(MutablePolynomial.createMonomial(1, 1)); // <- add x
        giantSteps.add(xPower);
        MutablePolynomial xPowerBig = xPower;
        int tBrentKung = (int) Math.sqrt(poly.degree);
        ArrayList<MutablePolynomial> hPowers = polyPowers(xPowerBig, poly, invMod, modulus, tBrentKung);
        for (int i = 0; i < m - 1; ++i)
            giantSteps.add(xPowerBig = compositionBrentKung(xPowerBig, hPowers, poly, invMod, modulus, tBrentKung));

        return new BabyGiantSteps(B, l, m, babySteps, giantSteps, invMod);
    }

    /** Shoup's main gcd loop */
    static FactorDecomposition DistinctDegreeFactorizationShoup(MutablePolynomial poly,
                                                                long modulus,
                                                                BabyGiantSteps steps) {
        //generate each I_j
        ArrayList<MutablePolynomial> iBases = new ArrayList<>();
        for (int j = 0; j <= steps.m; ++j) {
            MutablePolynomial iBase = MutablePolynomial.one();
            for (int i = 0; i <= steps.l - 1; ++i) {
                MutablePolynomial tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus);
                iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, modulus, false);
            }
            iBases.add(iBase);
        }

        ArrayList<MutablePolynomial> factors = new ArrayList<>();
        TIntArrayList degrees = new TIntArrayList();

        MutablePolynomial current = poly.clone();
        for (int j = 1; j <= steps.m; ++j) {
            MutablePolynomial gcd = PolynomialGCD(current, iBases.get(j), modulus);
            if (gcd.isConstant())
                continue;
            current = quotient(current, gcd, modulus, false);
            for (int i = steps.l - 1; i >= 0; --i) {
                MutablePolynomial tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus), modulus);
                if (!tmp.isConstant()) {
                    factors.add(tmp.clone().monic(modulus));
                    degrees.add(steps.l * j - i);
                }
                gcd = quotient(gcd, tmp, modulus, false);
            }
        }
        if (!current.isOne()) {
            factors.add(current.monic(modulus));
            degrees.add(current.degree);
        }

        return new FactorDecomposition(factors, degrees, 1);
    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} modulo {@code modulus} using
     * Victor Shoup's baby step / giant step algorithm.
     * <p>
     * In the case of not square-free input, the algorithm works, but the resulting d.d.f. may not be complete.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    public static FactorDecomposition DistinctDegreeFactorizationShoup(MutablePolynomial poly, long modulus) {
        long factor = mod(poly.lc(), modulus);
        poly = poly.clone().monic(modulus);
        return DistinctDegreeFactorizationShoup(poly, modulus, generateBabyGiantSteps(poly, modulus)).setFactor(factor);
    }

    /** baby/giant steps for Shoup's d.d.f. algorithm */
    static final class BabyGiantSteps {
        final int B, l, m;
        final ArrayList<MutablePolynomial> babySteps;
        final ArrayList<MutablePolynomial> giantSteps;
        final InverseModMonomial invMod;

        public BabyGiantSteps(int b, int l, int m, ArrayList<MutablePolynomial> babySteps, ArrayList<MutablePolynomial> giantSteps, InverseModMonomial invMod) {
            this.B = b;
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    public static FactorDecomposition DistinctDegreeFactorization(MutablePolynomial poly, long modulus) {
        if (poly.degree < 50)
            return DistinctDegreeFactorizationPlain(poly, modulus);
        else if (poly.degree < 200)
            return DistinctDegreeFactorizationPrecomputedExponents(poly, modulus);
        else
            return DistinctDegreeFactorizationShoup(poly, modulus);
    }

    static FactorDecomposition fixDistinctDegreeDecomposition(FactorDecomposition f) {
        for (int i = f.factors.size() - 1; i >= 0; --i) {
            int exponent = f.exponents.get(i);
            int degree = f.factors.get(i).degree;

            if (degree % exponent != 0)
                f.exponents.set(i, degree);
        }
        return f;
    }

    /**
     * Performs square-free factorization followed by distinct-degree factorization modulo {@code modulus}.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return square-free and distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition DistinctDegreeFactorizationComplete(MutablePolynomial poly, long modulus) {
        FactorDecomposition squareFree = SquareFreeFactorization(poly, modulus);
        long overallFactor = squareFree.factor;
        ArrayList<MutablePolynomial> finalFactors = new ArrayList<>(squareFree.factors.size());
        TIntArrayList finalExponents = new TIntArrayList(squareFree.factors.size());
        for (int i = squareFree.factors.size() - 1; i >= 0; --i) {
            FactorDecomposition dd = DistinctDegreeFactorization(squareFree.factors.get(i), modulus);
            int nFactors = dd.factors.size();
            finalFactors.ensureCapacity(finalFactors.size() + nFactors);
            finalExponents.ensureCapacity(finalExponents.size() + nFactors);
            for (int j = nFactors - 1; j >= 0; --j) {
                finalFactors.add(dd.factors.get(j));
                finalExponents.add(squareFree.exponents.get(i));
            }
            overallFactor = multiplyMod(overallFactor, dd.factor, modulus);
        }

        return new FactorDecomposition(finalFactors, finalExponents, overallFactor);
    }
//
//    static final class GiantSteps {
//        final ArrayList<MutablePolynomial> giantSteps = new ArrayList<>();
//        final MutablePolynomial poly;
//        final InverseModMonomial invMod;
//        final long modulus;
//        final MutablePolynomial basePower;
//
//        final int tBrentKung;
//        final ArrayList<MutablePolynomial> hPowers = new ArrayList<>();
//
//        public GiantSteps(MutablePolynomial poly, InverseModMonomial invMod, long modulus, MutablePolynomial basePower) {
//            this.poly = poly;
//            this.invMod = invMod;
//            this.modulus = modulus;
//            this.basePower = basePower;
//            this.tBrentKung = (int) Math.sqrt(poly.degree);
//
//            giantSteps.add(MutablePolynomial.createMonomial(1, 1)); // <- add x
//            giantSteps.add(basePower); // <- add x^p mod poly
//        }
//
//        MutablePolynomial get(int j) {
//            if (j < giantSteps.size())
//                return giantSteps.get(j);
//            if (hPowers.isEmpty())
//                polyPowers(basePower, poly, invMod, modulus, tBrentKung, hPowers);
//
//            MutablePolynomial xPowerBig = giantSteps.get(giantSteps.size() - 1);
//            for (int i = giantSteps.size(); i <= j; ++i)
//                giantSteps.add(xPowerBig = compositionBrentKung(xPowerBig, hPowers, poly, invMod, modulus, tBrentKung));
//            return xPowerBig;
//        }
//    }
//
//    static final class DDFSteps2 {
//        final int B, l, m;
//        final ArrayList<MutablePolynomial> babySteps;
//        final GiantSteps giantSteps;
//        final InverseModMonomial invMod;
//
//        public DDFSteps2(int b, int l, int m, ArrayList<MutablePolynomial> babySteps, GiantSteps giantSteps, InverseModMonomial invMod) {
//            B = b;
//            this.l = l;
//            this.m = m;
//            this.babySteps = babySteps;
//            this.giantSteps = giantSteps;
//            this.invMod = invMod;
//        }
//    }
//
//    static final class IBases {
//        final DDFSteps2 steps;
//        final long modulus;
//        final MutablePolynomial poly;
//
//        public IBases(DDFSteps2 steps, long modulus, MutablePolynomial poly) {
//            this.steps = steps;
//            this.modulus = modulus;
//            this.poly = poly;
//        }
//
//        final ArrayList<MutablePolynomial> iBases = new ArrayList<>();
//
//        public MutablePolynomial get(int k) {
//            if (k < iBases.size())
//                return iBases.get(k);
//
//            for (int j = iBases.size(); j <= k; ++j) {
//                MutablePolynomial iBase = MutablePolynomial.one();
//                for (int i = 0; i <= steps.l - 1; ++i) {
//                    MutablePolynomial tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus);
//                    iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, modulus, false);
//                }
//                iBases.add(iBase);
//            }
//            return iBases.get(k);
//        }
//    }
//
//    static DDFSteps2 ShoupDDF2(MutablePolynomial poly, long modulus) {
//        int n = poly.degree;
//        int B = (int) Math.floor(n / 2.);
//        int l = (int) Math.floor(Math.sqrt(B));
//        int m = (int) Math.ceil(1.0 * B / l);
//
//        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
//        ArrayList<MutablePolynomial> xPowers = xPowers(poly, invMod, modulus);
//
//        //baby steps
//        ArrayList<MutablePolynomial> babySteps = new ArrayList<>();
//        babySteps.add(MutablePolynomial.createMonomial(1, 1)); // <- add x
//        MutablePolynomial xPower = xPowers.get(1); // x^p mod poly
//        babySteps.add(xPower); // <- add x^p mod poly
//        for (int i = 0; i <= l - 2; ++i)
//            babySteps.add(xPower = powModulusMod(xPower, poly, invMod, modulus, xPowers));
//
//        // <- xPower = x^(p^l) mod poly
//        return new DDFSteps2(B, l, m, babySteps, new GiantSteps(poly, invMod, modulus, xPower), invMod);
//    }
//
//    static FactorDecomposition ShoupDDFFactorization2(MutablePolynomial poly, long modulus) {
//        ArrayList<MutablePolynomial> r = ShoupDDFFactorization2(poly, modulus, ShoupDDF2(poly, modulus));
//        return new FactorDecomposition(r, new TIntArrayList(), 1);
//    }
//
//
//    private static ArrayList<MutablePolynomial> ShoupDDFFactorization2(MutablePolynomial poly, long modulus, DDFSteps2 steps) {
//        IBases iBases = new IBases(steps, modulus, poly);
//
//        ArrayList<MutablePolynomial> fList = new ArrayList<>();
//        for (int i = 0; i <= poly.degree; ++i)
//            fList.add(MutablePolynomial.one());
//
//        MutablePolynomial current = poly.clone();
//        for (int j = 1; j <= steps.m; ++j) {
//            MutablePolynomial gcd = PolynomialGCD(current, iBases.get(j), modulus);
//            if (gcd.isConstant())
//                continue;
//            current = quotient(current, gcd, modulus, false);
//            for (int i = steps.l - 1; i >= 0; --i) {
//                MutablePolynomial tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus), modulus);
//                fList.set(steps.l * j - i, tmp);
//                gcd = quotient(gcd, tmp, modulus, false);
//            }
//        }
//        if (!current.isOne())
//            fList.set(current.degree - 1, current);
//
//        for (int i = fList.size() - 1; i >= 0; --i)
//            if (fList.get(i).isOne())
//                fList.remove(i);
//
//        return fList;
//    }
}
