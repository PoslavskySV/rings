package cc.r2.core.poly2;

import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.poly2.DivideAndRemainder.divideAndRemainder;
import static cc.r2.core.poly2.FactorDecomposition.oneFactor;
import static cc.r2.core.poly2.PolynomialGCD.PolynomialGCD;


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
    public static FactorDecomposition<MutablePolynomialZ> SquareFreeFactorization(MutablePolynomialZ poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun(poly);

        MutablePolynomialZ expFree = MutablePolynomialZ.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
        return SquareFreeFactorizationYun(expFree).addFactor(MutablePolynomialZ.createMonomial(1, exponent), 1);
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition<MutablePolynomialZ> SquareFreeFactorizationYun(MutablePolynomialZ poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutablePolynomialZ derivative = poly.derivative(), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutablePolynomialZ
                quot = divideAndRemainder(poly, gcd, false)[0], // safely destroy (cloned) poly (not used further)
                dQuot = divideAndRemainder(derivative, gcd, false)[0]; // safely destroy (cloned) derivative (not used further)

        ArrayList<MutablePolynomialZ> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            dQuot = dQuot.subtract(quot.derivative());
            MutablePolynomialZ factor = PolynomialGCD(quot, dQuot);
            quot = divideAndRemainder(quot, factor, false)[0]; // can destroy quot in divideAndRemainder
            dQuot = divideAndRemainder(dQuot, factor, false)[0]; // can destroy dQuot in divideAndRemainder
            if (!factor.isOne()) {
                factors.add(factor);
                exponents.add(i);
            }
        }
        return new FactorDecomposition<>(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a poly using Musser's algorithm
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static FactorDecomposition<MutablePolynomialZ> SquareFreeFactorizationMusser(MutablePolynomialZ poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutablePolynomialZ derivative = poly.derivative(), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutablePolynomialZ quot = divideAndRemainder(poly, gcd, false)[0]; // safely destroy (cloned) poly

        ArrayList<MutablePolynomialZ> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (true) {
            ++i;
            MutablePolynomialZ nextQuot = PolynomialGCD(gcd, quot);
            gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // safely destroy gcd (reassigned)
            MutablePolynomialZ factor = divideAndRemainder(quot, nextQuot, false)[0]; // safely destroy quot (reassigned further)
            if (!factor.isConstant()) {
                factors.add(factor);
                exponents.add(i);
            }
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
        return new FactorDecomposition<>(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a {@code poly} modulo {@code modulus}.
     *
     * @param poly the polynomial
     * @return square-free decomposition modulo {@code modulus}
     */
    public static FactorDecomposition<MutablePolynomialMod> SquareFreeFactorization(MutablePolynomialMod poly) {
//        if (modulus >= Integer.MAX_VALUE)// <- just to be on the safe side)
//            throw new IllegalArgumentException();

        poly = poly.clone();
        long lc = poly.lc();
        //make poly monic
        poly = poly.monic();

        if (poly.isConstant())
            return oneFactor(lc);

        if (poly.degree <= 1)
            return oneFactor(poly, lc);

        FactorDecomposition<MutablePolynomialMod> factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly);
        else {
            MutablePolynomialMod expFree = MutablePolynomialMod.create(poly.modulus, Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
            factorization = SquareFreeFactorizationMusser0(expFree).addFactor(MutablePolynomialMod.createMonomial(poly.modulus, 1, exponent), 1);
        }

        return factorization.setFactor(lc);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static FactorDecomposition<MutablePolynomialMod> SquareFreeFactorizationMusser0(MutablePolynomialMod poly) {
        poly.monic();
        if (poly.isConstant())
            return oneFactor(poly.lc());

        if (poly.degree <= 1)
            return oneFactor(poly, 1);

        MutablePolynomialMod derivative = poly.derivative();
        if (!derivative.isZero()) {
            MutablePolynomialMod gcd = PolynomialGCD(poly, derivative);
            if (gcd.isConstant())
                return oneFactor(poly, 1);
            MutablePolynomialMod quot = divideAndRemainder(poly, gcd, false)[0]; // can safely destroy poly (not used further)

            ArrayList<MutablePolynomialMod> factors = new ArrayList<>();
            TIntArrayList exponents = new TIntArrayList();
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                MutablePolynomialMod nextQuot = PolynomialGCD(gcd, quot);
                MutablePolynomialMod factor = divideAndRemainder(quot, nextQuot, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant()) {
                    factors.add(factor.monic());
                    exponents.add(i);
                }
                gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd);
                FactorDecomposition<MutablePolynomialMod> gcdFactorization =
                        SquareFreeFactorizationMusser0(gcd)
                                .raiseExponents(poly.modulus);

                factors.addAll(gcdFactorization.factors);
                exponents.addAll(gcdFactorization.exponents);

                return new FactorDecomposition<>(factors, exponents, 1);
            } else
                return new FactorDecomposition<>(factors, exponents, 1);
        } else {
            MutablePolynomialMod pRoot = pRoot(poly);
            return SquareFreeFactorizationMusser0(pRoot).raiseExponents(poly.modulus).setFactor(1);
        }
    }

    /** p-th root of poly */
    private static MutablePolynomialMod pRoot(MutablePolynomialMod poly) {
        long modulus = poly.modulus;
        assert poly.degree % modulus == 0;
        long[] rootData = new long[poly.degree / (int) modulus + 1];
        Arrays.fill(rootData, 0);
        for (int i = poly.degree; i >= 0; --i)
            if (poly.data[i] != 0) {
                assert i % modulus == 0;
                rootData[i / (int) modulus] = poly.data[i];
            }
        return poly.createFromArray(rootData);
    }

    /**
     * Returns {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     */
    public static <T extends MutablePolynomialAbstract<T>> boolean isSquareFree(T poly) {
        return PolynomialGCD(poly, poly.derivative()).isConstant();
    }

    @SuppressWarnings("unchecked")
    public static <T extends MutablePolynomialAbstract<T>> FactorDecomposition<T> SquareFreeFactorization(T poly) {
        if (poly instanceof MutablePolynomialZ)
            return (FactorDecomposition<T>) SquareFreeFactorization((MutablePolynomialZ) poly);
        else
            return (FactorDecomposition<T>) SquareFreeFactorization((MutablePolynomialMod) poly);
    }

    /**
     * Returns square-free part of the {@code poly}
     *
     * @param poly the polynomial
     * @return square free part
     */
    public static <T extends MutablePolynomialAbstract<T>> T SquareFreePart(T poly) {
        return SquareFreeFactorization(poly).factors.stream().filter(x -> !x.isMonomial()).reduce(poly.createOne(), MutablePolynomialAbstract<T>::multiply);
    }
}