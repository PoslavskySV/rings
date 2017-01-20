package cc.r2.core.polynomial;

import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;

import static cc.r2.core.polynomial.DivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.DivideAndRemainder.quotient;
import static cc.r2.core.polynomial.Factorization.oneFactor;
import static cc.r2.core.polynomial.LongArithmetics.*;
import static cc.r2.core.polynomial.ModularComposition.*;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyMultiplyMod;
import static cc.r2.core.polynomial.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.polynomial.SquareFreeFactorization.SquareFreeFactorization;

/**
 * Created by poslavsky on 20/01/2017.
 */
public final class DistinctDegreeFactorization {
    private DistinctDegreeFactorization() {}

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} modulo {@code modulus}.
     * In the case of not square-free input, the algorithm works, but the resulting factorization is not compete.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    static Factorization DistinctDegreeFactorizationPlain(
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

        DivideAndRemainder.InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
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
                    degrees.add(i + 1);
                }
                break;
            }
        }

        return new Factorization(factors, degrees, factor);
    }

    static BabyGiantSteps generateBabyGiantSteps(MutablePolynomial poly, long modulus) {
        int n = poly.degree;
        int B = (int) Math.floor(n / 2.);
        int l = (int) Math.floor(Math.sqrt(B));
        int m = (int) Math.ceil(1.0 * B / l);

        DivideAndRemainder.InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
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

    static Factorization DistinctDegreeFactorizationShoup(MutablePolynomial poly,
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
                factors.add(tmp);
                degrees.add(steps.l * j - i);
                gcd = quotient(gcd, tmp, modulus, false);
            }
        }
        if (!current.isOne()) {
            factors.add(current);
            degrees.add(current.degree - 1);
        }

        return new Factorization(factors, degrees, 1);
    }

    static Factorization DistinctDegreeFactorizationShoup(MutablePolynomial poly,
                                                          long modulus) {
        return DistinctDegreeFactorizationShoup(poly, modulus, generateBabyGiantSteps(poly, modulus));
    }

    /** baby-giant steps for Shoup's algorithm for DDF */
    static final class BabyGiantSteps {
        final int B, l, m;
        final ArrayList<MutablePolynomial> babySteps;
        final ArrayList<MutablePolynomial> giantSteps;
        final DivideAndRemainder.InverseModMonomial invMod;

        public BabyGiantSteps(int b, int l, int m, ArrayList<MutablePolynomial> babySteps, ArrayList<MutablePolynomial> giantSteps, DivideAndRemainder.InverseModMonomial invMod) {
            B = b;
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    /**
     * Performs square-free factorization followed by distinct-degree factorization modulo {@code modulus}.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return square-free and distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    public static Factorization DistinctDegreeFactorizationComplete(MutablePolynomial poly, long modulus) {
        Factorization squareFree = SquareFreeFactorization(poly, modulus);
        long overallFactor = squareFree.factor;
        ArrayList<MutablePolynomial> finalFactors = new ArrayList<>(squareFree.factors.length);
        TIntArrayList finalExponents = new TIntArrayList(squareFree.factors.length);
        for (int i = squareFree.factors.length - 1; i >= 0; --i) {
            Factorization dd = DistinctDegreeFactorizationPlain(squareFree.factors[i], modulus);
            int nFactors = dd.factors.length;
            finalFactors.ensureCapacity(finalFactors.size() + nFactors);
            finalExponents.ensureCapacity(finalExponents.size() + nFactors);
            for (int j = nFactors - 1; j >= 0; --j) {
                finalFactors.add(dd.factors[j]);
                finalExponents.add(squareFree.exponents[i]);
            }
            overallFactor = multiplyMod(overallFactor, dd.factor, modulus);
        }

        return new Factorization(finalFactors, finalExponents, overallFactor);
    }
}
