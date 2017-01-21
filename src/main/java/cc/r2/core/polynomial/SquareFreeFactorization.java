package cc.r2.core.polynomial;

import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.DivideAndRemainder.divideAndRemainder;
import static cc.r2.core.polynomial.FactorDecomposition.oneFactor;
import static cc.r2.core.polynomial.MutablePolynomial.derivative;
import static cc.r2.core.polynomial.PolynomialGCD.PolynomialGCD;

/**
 * Square-free factorization of univariate polynomials with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 */
public final class SquareFreeFactorization {
    private SquareFreeFactorization() {}

    /**
     * Performs square-free factorization of a poly.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static FactorDecomposition SquareFreeFactorization(MutablePolynomial poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun(poly);

        MutablePolynomial expFree = MutablePolynomial.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
        return SquareFreeFactorizationYun(expFree).addFactor(MutablePolynomial.createMonomial(1, exponent), 1);
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition SquareFreeFactorizationYun(MutablePolynomial poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divide(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutablePolynomial derivative = derivative(poly), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutablePolynomial
                quot = divideAndRemainder(poly, gcd, false)[0], // safely destroy (cloned) poly (not used further)
                dQuot = divideAndRemainder(derivative, gcd, false)[0]; // safely destroy (cloned) derivative (not used further)

        ArrayList<MutablePolynomial> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            dQuot = dQuot.subtract(derivative(quot));
            MutablePolynomial factor = PolynomialGCD(quot, dQuot);
            quot = divideAndRemainder(quot, factor, false)[0]; // can destroy quot in divideAndRemainder
            dQuot = divideAndRemainder(dQuot, factor, false)[0]; // can destroy dQuot in divideAndRemainder
            if (!factor.isOne()) {
                factors.add(factor);
                exponents.add(i);
            }
        }
        return new FactorDecomposition(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a poly using Musser's algorithm
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition SquareFreeFactorizationMusser(MutablePolynomial poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divide(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutablePolynomial derivative = derivative(poly), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutablePolynomial quot = divideAndRemainder(poly, gcd, false)[0]; // safely destroy (cloned) poly

        ArrayList<MutablePolynomial> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (true) {
            ++i;
            MutablePolynomial nextQuot = PolynomialGCD(gcd, quot);
            gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // safely destroy gcd (reassigned)
            MutablePolynomial factor = divideAndRemainder(quot, nextQuot, false)[0]; // safely destroy quot (reassigned further)
            if (!factor.isConstant()) {
                factors.add(factor);
                exponents.add(i);
            }
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
        return new FactorDecomposition(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a {@code poly} modulo {@code modulus}.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return square-free decomposition modulo {@code modulus}
     */
    public static FactorDecomposition SquareFreeFactorization(MutablePolynomial poly, long modulus) {
        if (modulus >= Integer.MAX_VALUE)// <- just to be on the safe side)
            throw new IllegalArgumentException();

        poly = poly.clone().modulus(modulus);
        long lc = poly.lc();
        //make poly monic
        poly = poly.monic(modulus);

        if (poly.isConstant())
            return oneFactor(lc);

        if (poly.degree <= 1)
            return oneFactor(poly, lc);

        FactorDecomposition factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly, modulus);
        else {
            MutablePolynomial expFree = MutablePolynomial.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
            factorization = SquareFreeFactorizationMusser0(expFree, modulus).addFactor(MutablePolynomial.createMonomial(1, exponent), 1);
        }

        return factorization.setFactor(lc);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static FactorDecomposition SquareFreeFactorizationMusser0(MutablePolynomial poly, long modulus) {
        poly.monic(modulus);
        if (poly.isConstant())
            return oneFactor(poly.lc());

        if (poly.degree <= 1)
            return oneFactor(poly, 1);

        MutablePolynomial derivative = derivative(poly, modulus);
        if (!derivative.isZero()) {
            MutablePolynomial gcd = PolynomialGCD(poly, derivative, modulus);
            if (gcd.isConstant())
                return oneFactor(poly, 1);
            MutablePolynomial quot = divideAndRemainder(poly, gcd, modulus, false)[0]; // can safely destroy poly (not used further)

            ArrayList<MutablePolynomial> factors = new ArrayList<>();
            TIntArrayList exponents = new TIntArrayList();
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                MutablePolynomial nextQuot = PolynomialGCD(gcd, quot, modulus);
                MutablePolynomial factor = divideAndRemainder(quot, nextQuot, modulus, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant()) {
                    factors.add(factor.monic(modulus));
                    exponents.add(i);
                }
                gcd = divideAndRemainder(gcd, nextQuot, modulus, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd, modulus);
                FactorDecomposition gcdFactorization =
                        SquareFreeFactorizationMusser0(gcd, modulus)
                                .raiseExponents(modulus);

                factors.addAll(gcdFactorization.factors);
                exponents.addAll(gcdFactorization.exponents);

                return new FactorDecomposition(factors, exponents, 1);
            } else
                return new FactorDecomposition(factors, exponents, 1);
        } else {
            MutablePolynomial pRoot = pRoot(poly, modulus);
            return SquareFreeFactorizationMusser0(pRoot, modulus).raiseExponents(modulus).setFactor(1);
        }
    }

    /** p-th root of poly */
    private static MutablePolynomial pRoot(MutablePolynomial poly, long modulus) {
        assert poly.degree % modulus == 0;
        long[] rootData = new long[poly.degree / (int) modulus + 1];
        Arrays.fill(rootData, 0);
        for (int i = poly.degree; i >= 0; --i)
            if (poly.data[i] != 0) {
                assert i % modulus == 0;
                rootData[i / (int) modulus] = poly.data[i];
            }
        return MutablePolynomial.create(rootData);
    }

    /**
     * Returns {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     */
    public static boolean isSquareFree(MutablePolynomial poly) {
        return PolynomialGCD(poly, derivative(poly)).isConstant();
    }

    /**
     * Returns {@code true} if and only if {@code poly} is square-free modulo {@code modulus} and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if and only if {@code poly} is square-free modulo {@code modulus} and {@code false} otherwise
     */
    public static boolean isSquareFree(MutablePolynomial poly, long modulus) {
        return PolynomialGCD(poly, derivative(poly, modulus), modulus).isConstant();
    }

    /**
     * Returns square-free part of the {@code poly}
     *
     * @param poly the polynomial
     * @return square free part
     */
    public static MutablePolynomial SquareFreePart(MutablePolynomial poly) {
        return SquareFreeFactorization(poly).factors.stream().filter(x -> !x.isMonomial()).reduce(MutablePolynomial.one(), MutablePolynomial::multiply);
    }

    /**
     * Returns square-free part of the {@code poly} modulo {@code modulus}
     *
     * @param poly    the polynomial
     * @param modulus the modulus
     * @return square free part modulo {@code modulus}
     */
    public static MutablePolynomial SquareFreePart(MutablePolynomial poly, long modulus) {
        return SquareFreeFactorization(poly, modulus).factors.stream().filter(x -> !x.isMonomial()).reduce(MutablePolynomial.one(), (a, b) -> a.multiply(b, modulus));
    }
}