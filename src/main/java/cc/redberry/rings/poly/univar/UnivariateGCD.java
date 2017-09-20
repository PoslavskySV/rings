package cc.redberry.rings.poly.univar;


import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.bigint.ChineseRemainders;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.primes.PrimesIterator;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;

import static cc.redberry.rings.bigint.ChineseRemainders.ChineseRemainders;
import static cc.redberry.rings.poly.univar.Conversions64bit.asOverZp64;
import static cc.redberry.rings.poly.univar.Conversions64bit.canConvertToZp64;

/**
 * Univariate polynomial GCD and sub-resultant sequences.
 *
 * @since 1.0
 */
public final class UnivariateGCD {
    private UnivariateGCD() {}

    /**
     * Calculates the GCD of two polynomials. Depending on the coefficient ring, the algorithm switches between Half-GCD
     * (polys over finite fields), modular GCD (polys over Z and Q) and subresultant Euclid (other rings).
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> T PolynomialGCD(T a, T b) {
        a.assertSameCoefficientRingWith(b);
        if (a.isOverField())
            return HalfGCD(a, b);
        else if (a instanceof UnivariatePolynomialZ64)
            return (T) ModularGCD((UnivariatePolynomialZ64) a, (UnivariatePolynomialZ64) b);
        else if (a instanceof UnivariatePolynomial) {
            Ring ring = ((UnivariatePolynomial) a).ring;
            if (ring.equals(Rings.Z))
                return (T) ModularGCD((UnivariatePolynomial) a, (UnivariatePolynomial) b);
            else
                return (T) EuclidSubresultantRemainders((UnivariatePolynomial) a, (UnivariatePolynomial) b).gcd();
        } else
            throw new RuntimeException(a.getClass().toString());
    }

    /**
     * Computes {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}. Half-GCD algorithm is used.
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)} (gcd is monic)
     * @see #ExtendedHalfGCD(IUnivariatePolynomial, IUnivariatePolynomial)
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> T[] PolynomialExtendedGCD(T a, T b) {
        if (a.isOverField())
            return ExtendedHalfGCD(a, b);
        else
            throw new IllegalArgumentException("Polynomial over field is expected");
    }

    /**
     * Returns GCD of a list of polynomials.
     *
     * @param polynomials a set of polynomials
     * @return GCD of polynomials
     */
    public static <T extends IUnivariatePolynomial<T>> T PolynomialGCD(T... polynomials) {
        T gcd = polynomials[0];
        for (int i = 1; i < polynomials.length; i++)
            gcd = PolynomialGCD(gcd, polynomials[i]);
        return gcd;
    }

    /**
     * Returns GCD of a list of polynomials.
     *
     * @param polynomials a set of polynomials
     * @return GCD of polynomials
     */
    public static <T extends IUnivariatePolynomial<T>> T PolynomialGCD(Iterable<T> polynomials) {
        T gcd = null;
        for (T poly : polynomials)
            gcd = gcd == null ? poly : PolynomialGCD(gcd, poly);
        return gcd;
    }

    /* ========================================== implementation ==================================================== */

    @SuppressWarnings("unchecked")
    static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> convert(PolynomialRemainders<UnivariatePolynomialZp64> r64) {
        return new PolynomialRemainders<>(new ArrayList<>(r64.remainders.stream().map(p -> (T) p.toBigPoly()).collect(Collectors.toList())));
    }

    /**
     * Returns the remainder sequence produced by Euclidean algorithm. If coefficient ring of the input is not a field
     * (and thus polynomials does not form an integral ring), {@code ArithmeticException} may be thrown in case when
     * some exact divisions are not possible.
     *
     * @param a poly
     * @param b poly
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> EuclidRemainders(final T a, final T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            // switch to faster machine arithmetic
            return convert(EuclidRemainders(asOverZp64(a), asOverZp64(b)));

        if (a.degree() < b.degree())
            return EuclidRemainders(b, a);

        ArrayList<T> prs = new ArrayList<>();
        prs.add(a.clone());
        prs.add(b.clone());

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(prs);

        T x = a, y = b, r;
        while (true) {
            r = UnivariateDivision.remainder(x, y, true);
            if (r == null)
                throw new ArithmeticException("Not divisible with remainder: (" + x + ") / (" + y + ")");

            if (r.isZero())
                break;
            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders<>(prs);
    }

    /**
     * Returns the GCD calculated with Euclidean algorithm. If coefficient ring of the input is not a field (and thus
     * polynomials does not form an integral ring), {@code ArithmeticException} may be thrown in case when some exact
     * divisions are not possible.
     *
     * @param a poly
     * @param b poly
     * @return the GCD (monic if a and b are over field)
     */
    public static <T extends IUnivariatePolynomial<T>> T EuclidGCD(final T a, final T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            return Conversions64bit.convert(EuclidGCD(asOverZp64(a), asOverZp64(b)));

        if (a.degree() < b.degree())
            return EuclidGCD(b, a);

        if (a.isZero()) return normalizeGCD(b.clone());
        if (b.isZero()) return normalizeGCD(a.clone());

        T x = a, y = b, r;
        while (true) {
            r = UnivariateDivision.remainder(x, y, true);
            if (r == null)
                throw new ArithmeticException("Not divisible with remainder: (" + x + ") / (" + y + ")");

            if (r.isZero())
                break;
            x = y;
            y = r;
        }
        return normalizeGCD(y == a ? y.clone() : (y == b ? y.clone() : y));
    }

    private static <T extends IUnivariatePolynomial<T>> T normalizeGCD(T gcd) {
        if (gcd.isOverField())
            return gcd.monic();
        else
            return gcd;
    }

    /**
     * Runs extended Euclidean algorithm to compute {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a,
     * b)}. If coefficient ring of the input is not a field (and thus polynomials does not form an integral ring),
     * {@code ArithmeticException} may be thrown in case when some exact divisions are not possible.
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> T[] ExtendedEuclidGCD(final T a, final T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            return Conversions64bit.convert(a, ExtendedEuclidGCD(asOverZp64(a), asOverZp64(b)));

        T s = a.createZero(), old_s = a.createOne();
        T t = a.createOne(), old_t = a.createZero();
        T r = b.clone(), old_r = a.clone();

        T q;
        T tmp;
        while (!r.isZero()) {
            q = UnivariateDivision.quotient(old_r, r, true);
            if (q == null)
                throw new ArithmeticException("Not divisible with remainder: (" + old_r + ") / (" + r + ")");

            tmp = old_r;
            old_r = r;
            r = tmp.clone().subtract(q.clone().multiply(r));

            tmp = old_s;
            old_s = s;
            s = tmp.clone().subtract(q.clone().multiply(s));

            tmp = old_t;
            old_t = t;
            t = tmp.clone().subtract(q.clone().multiply(t));
        }
        assert old_r.equals(a.clone().multiply(old_s).add(b.clone().multiply(old_t))) : a.clone().multiply(old_s).add
                (b.clone().multiply(old_t));

        T[] result = a.createArray(3);
        result[0] = old_r;
        result[1] = old_s;
        result[2] = old_t;
        return normalizeExtendedGCD(result);
    }

    private static <T extends IUnivariatePolynomial<T>> T[] normalizeExtendedGCD(T[] xgcd) {
        if (!xgcd[0].isOverField())
            return xgcd;

        if (xgcd[0].isZero())
            return xgcd;
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[0].monic();
        return xgcd;
    }

    /** for polynomial degrees larger than this a Half-GCD algorithm will be used */
    static int SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE = 180;
    /** for polynomial degrees larger than this a Half-GCD algorithm for hMatrix will be used */
    static int SWITCH_TO_HALF_GCD_H_MATRIX_DEGREE = 25;

    /**
     * Half-GCD algorithm. The algorithm automatically switches to Euclidean algorithm for small input. If coefficient
     * ring of the input is not a field (and thus polynomials does not form an integral ring), {@code
     * ArithmeticException} may be thrown in case when some exact divisions are not possible.
     *
     * @param a poly
     * @param b poly
     * @return the GCD (monic)
     */
    public static <T extends IUnivariatePolynomial<T>> T HalfGCD(T a, T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            return Conversions64bit.convert(HalfGCD(asOverZp64(a), asOverZp64(b)));

        if (a.degree() < b.degree())
            return HalfGCD(b, a);

        if (a.degree() == b.degree())
            b = UnivariateDivision.remainder(b, a, true);

        while (a.degree() > SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE && !b.isZero()) {
            T[] col = reduceHalfGCD(a, b);
            a = col[0];
            b = col[1];

            if (!b.isZero()) {
                T remainder = UnivariateDivision.remainder(a, b, true);
                if (remainder == null)
                    throw new ArithmeticException("Not divisible with remainder: (" + a + ") / (" + b + ")");
                a = b;
                b = remainder;
            }
        }
        return UnivariateGCD.EuclidGCD(a, b);
    }

    /**
     * Runs extended Half-GCD algorithm to compute {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}.
     * If coefficient ring of the input is not a field (and thus polynomials does not form an integral ring), {@code
     * ArithmeticException} may be thrown in case when some exact divisions are not possible.
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)} (gcd is monic)
     */
    public static <T extends IUnivariatePolynomial<T>> T[] ExtendedHalfGCD(T a, T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            return Conversions64bit.convert(a, ExtendedHalfGCD(asOverZp64(a), asOverZp64(b)));

        if (a.degree() < b.degree()) {
            T[] r = ExtendedHalfGCD(b, a);
            ArraysUtil.swap(r, 1, 2);
            return r;
        }

        if (b.isZero()) {
            T[] result = a.createArray(3);
            result[0] = a.clone();
            result[1] = a.createOne();
            result[2] = a.createZero();
            return normalizeExtendedGCD(result);
        }

        a = a.clone();
        b = b.clone();

        T quotient = null;
        if (a.degree() == b.degree()) {
            T[] qd = UnivariateDivision.divideAndRemainder(a, b, true);
            if (qd == null)
                throw new ArithmeticException("Not divisible with remainder: (" + a + ") / (" + b + ")");
            quotient = qd[0];
            T remainder = qd[1];
            a = b;
            b = remainder;
        }

        T[][] hMatrix = reduceExtendedHalfGCD(a, b, a.degree() + 1);
        T gcd = a, s, t;

        if (quotient != null) {
            s = hMatrix[0][1];
            t = quotient.multiply(hMatrix[0][1]);
            t = hMatrix[0][0].subtract(t);
        } else {
            s = hMatrix[0][0];
            t = hMatrix[0][1];
        }

        T[] result = a.createArray(3);
        result[0] = gcd;
        result[1] = s;
        result[2] = t;
        return normalizeExtendedGCD(result);
    }

    /**
     * @param reduce whether to reduce a and b inplace
     */
    private static <T extends IUnivariatePolynomial<T>> T[][] hMatrixPlain(T a, T b, int degreeToReduce, boolean reduce) {
        T[][] hMatrix = unitMatrix(a);
        int goal = a.degree() - degreeToReduce;
        if (b.degree() <= goal)
            return hMatrix;

        T tmpA = a, tmpB = b;
        while (tmpB.degree() > goal && !tmpB.isZero()) {
            T[] qd = UnivariateDivision.divideAndRemainder(tmpA, tmpB, true);
            if (qd == null)
                throw new ArithmeticException("Not divisible with remainder: (" + tmpA + ") / (" + tmpB + ")");
            T quotient = qd[0], remainder = qd[1];

            T tmp;
            tmp = quotient.clone().multiply(hMatrix[1][0]);
            tmp = hMatrix[0][0].clone().subtract(tmp);

            hMatrix[0][0] = hMatrix[1][0];
            hMatrix[1][0] = tmp;

            tmp = quotient.clone().multiply(hMatrix[1][1]);
            tmp = hMatrix[0][1].clone().subtract(tmp);

            hMatrix[0][1] = hMatrix[1][1];
            hMatrix[1][1] = tmp;

            tmpA = tmpB;
            tmpB = remainder;
        }

        if (reduce) {
            a.setAndDestroy(tmpA);
            b.setAndDestroy(tmpB);
        }

        return hMatrix;
    }

    private static <T extends IUnivariatePolynomial<T>> T[][] hMatrixHalfGCD(T a, T b, int d) {
        if (b.isZero() || b.degree() <= a.degree() - d)
            return unitMatrix(a);

        int n = a.degree() - 2 * d + 2;
        if (n < 0) n = 0;

        T a1 = a.clone().shiftLeft(n);
        T b1 = b.clone().shiftLeft(n);

        if (d <= SWITCH_TO_HALF_GCD_H_MATRIX_DEGREE)
            return hMatrixPlain(a1, b1, d, false);

        int dR = (d + 1) / 2;
        if (dR < 1)
            dR = 1;
        if (dR >= d)
            dR = d - 1;

        T[][] hMatrixR = hMatrixHalfGCD(a1, b1, dR);
        T[] col = columnMultiply(hMatrixR, a1, b1);
        a1 = col[0];
        b1 = col[1];


        int dL = b1.degree() - a.degree() + n + d;
        if (b1.isZero() || dL <= 0)
            return hMatrixR;

        T[] qd = UnivariateDivision.divideAndRemainder(a1, b1, false);
        if (qd == null)
            throw new ArithmeticException("Not divisible with remainder: (" + a1 + ") / (" + b1 + ")");
        T quotient = qd[0], remainder = qd[1];
        T[][] hMatrixL = hMatrixHalfGCD(b1, remainder, dL);

        T tmp;
        tmp = quotient.clone().multiply(hMatrixR[1][0]);
        tmp = hMatrixR[0][0].clone().subtract(tmp);
        hMatrixR[0][0] = hMatrixR[1][0];
        hMatrixR[1][0] = tmp;

        tmp = quotient.clone().multiply(hMatrixR[1][1]);
        tmp = hMatrixR[0][1].clone().subtract(tmp);
        hMatrixR[0][1] = hMatrixR[1][1];
        hMatrixR[1][1] = tmp;

        return matrixMultiply(hMatrixL, hMatrixR);
    }

    /** a and b will be modified */
    static <T extends IUnivariatePolynomial<T>> T[][] reduceExtendedHalfGCD(T a, T b, int d) {
        assert a.degree() >= b.degree();
        if (b.isZero() || b.degree() <= a.degree() - d)
            return unitMatrix(a);

        int aDegree = a.degree();
        if (d <= SWITCH_TO_HALF_GCD_H_MATRIX_DEGREE)
            return hMatrixPlain(a, b, d, true);

        int dL = (d + 1) / 2;
        if (dL < 1)
            dL = 1;
        if (dL >= d)
            dL = d - 1;

        T[][] hMatrixR = hMatrixHalfGCD(a, b, dL);
        T[] col = columnMultiply(hMatrixR, a, b);
        a.setAndDestroy(col[0]);
        b.setAndDestroy(col[1]);

        int dR = b.degree() - aDegree + d;
        if (b.isZero() || dR <= 0)
            return hMatrixR;

        T[] qd = UnivariateDivision.divideAndRemainder(a, b, true);
        if (qd == null)
            throw new ArithmeticException("Not divisible with remainder: (" + a + ") / (" + b + ")");
        T quotient = qd[0], remainder = qd[1];

        a.setAndDestroy(b);
        b.setAndDestroy(remainder);
        T[][] hMatrixL = reduceExtendedHalfGCD(a, b, dR);

        T tmp;
        tmp = quotient.clone().multiply(hMatrixR[1][0]);
        tmp = hMatrixR[0][0].clone().subtract(tmp);

        hMatrixR[0][0] = hMatrixR[1][0];
        hMatrixR[1][0] = tmp;


        tmp = quotient.clone().multiply(hMatrixR[1][1]);
        tmp = hMatrixR[0][1].clone().subtract(tmp);
        hMatrixR[0][1] = hMatrixR[1][1];
        hMatrixR[1][1] = tmp;

        return matrixMultiply(hMatrixL, hMatrixR);
    }

    /** a and b will be modified */
    static <T extends IUnivariatePolynomial<T>> T[] reduceHalfGCD(T a, T b) {
        int d = (a.degree() + 1) / 2;

        if (b.isZero() || b.degree() <= a.degree() - d)
            return a.createArray(a, b);

        int aDegree = a.degree();

        int d1 = (d + 1) / 2;
        if (d1 < 1)
            d1 = 1;
        if (d1 >= d)
            d1 = d - 1;


        T[][] hMatrix = hMatrixHalfGCD(a, b, d1);
        T[] col = columnMultiply(hMatrix, a, b);
        a = col[0];
        b = col[1];

        int d2 = b.degree() - aDegree + d;

        if (b.isZero() || d2 <= 0)
            return a.createArray(a, b);

        T remainder = UnivariateDivision.remainder(a, b, true);
        if (remainder == null)
            throw new ArithmeticException("Not divisible with remainder: (" + a + ") / (" + b + ")");
        a = b;
        b = remainder;

        return columnMultiply(hMatrixHalfGCD(a, b, d2), a, b);
    }

    private static <T extends IUnivariatePolynomial<T>> T[][] matrixMultiply(T[][] matrix1, T[][] matrix2) {
        T[][] r = matrix1[0][0].createArray2d(2, 2);
        r[0][0] = matrix1[0][0].clone().multiply(matrix2[0][0]).add(matrix1[0][1].clone().multiply(matrix2[1][0]));
        r[0][1] = matrix1[0][0].clone().multiply(matrix2[0][1]).add(matrix1[0][1].clone().multiply(matrix2[1][1]));
        r[1][0] = matrix1[1][0].clone().multiply(matrix2[0][0]).add(matrix1[1][1].clone().multiply(matrix2[1][0]));
        r[1][1] = matrix1[1][0].clone().multiply(matrix2[0][1]).add(matrix1[1][1].clone().multiply(matrix2[1][1]));
        return r;
    }

    private static <T extends IUnivariatePolynomial<T>> T[] columnMultiply(T[][] hMatrix, T row1, T row2) {
        T[] resultColumn = row1.createArray(2);
        resultColumn[0] = hMatrix[0][0].clone().multiply(row1).add(hMatrix[0][1].clone().multiply(row2));
        resultColumn[1] = hMatrix[1][0].clone().multiply(row1).add(hMatrix[1][1].clone().multiply(row2));
        return resultColumn;
    }

    private static <T extends IUnivariatePolynomial<T>> T[][] unitMatrix(T factory) {
        T[][] m = factory.createArray2d(2, 2);
        m[0][0] = factory.createOne();
        m[0][1] = factory.createZero();
        m[1][0] = factory.createZero();
        m[1][1] = factory.createOne();
        return m;
    }

    /**
     * Euclidean algorithm for polynomials over Z that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to build primitive polynomial remainders or not
     * @return polynomial remainder sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> EuclidPseudoRemainders(final T a,
                                                                                                      final T b,
                                                                                                      boolean primitivePRS) {
        if (a instanceof UnivariatePolynomialZ64)
            return (PolynomialRemainders<T>) EuclidPseudoRemainders((UnivariatePolynomialZ64) a, (UnivariatePolynomialZ64) b, primitivePRS);
        else if (a instanceof UnivariatePolynomial)
            return (PolynomialRemainders<T>) EuclidPseudoRemainders((UnivariatePolynomial) a, (UnivariatePolynomial) b, primitivePRS);
        else
            throw new RuntimeException("Not a Z[x] polynomials: " + a.getClass());
    }

    /**
     * Euclidean algorithm for polynomials over Euclidean rings that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to build primitive polynomial remainders or not
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static PolynomialRemainders<UnivariatePolynomialZ64> EuclidPseudoRemainders(final UnivariatePolynomialZ64 a,
                                                                                       final UnivariatePolynomialZ64 b,
                                                                                       boolean primitivePRS) {
        if (a.degree < b.degree)
            return EuclidPseudoRemainders(b, a, primitivePRS);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());

        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        UnivariatePolynomialZ64 aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);
        PolynomialRemainders<UnivariatePolynomialZ64> res = EuclidPseudoRemainders0(aPP, bPP, primitivePRS);
        res.gcd().primitivePartSameSign().multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials over Euclidean rings that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to build primitive polynomial remainders or not
     * @return polynomial remainder sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <E> PolynomialRemainders<UnivariatePolynomial<E>> EuclidPseudoRemainders(final UnivariatePolynomial<E> a,
                                                                                           final UnivariatePolynomial<E> b,
                                                                                           boolean primitivePRS) {
        if (a.degree < b.degree)
            return EuclidPseudoRemainders(b, a, primitivePRS);

        if (a.isZero() || b.isZero())
            return new PolynomialRemainders<>(a.clone(), b.clone());

        E aContent = a.content(), bContent = b.content();
        E contentGCD = a.ring.gcd(aContent, bContent);
        UnivariatePolynomial<E> aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);
        PolynomialRemainders<UnivariatePolynomial<E>> res = EuclidPseudoRemainders0(aPP, bPP, primitivePRS);
        res.gcd().primitivePartSameSign().multiply(contentGCD);
        return res;
    }

    @SuppressWarnings("unchecked")
    private static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> EuclidPseudoRemainders0(final T aPP,
                                                                                                        final T bPP,
                                                                                                        boolean primitivePRS) {
        ArrayList<T> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        T x = aPP, y = bPP, r;
        while (true) {
            T[] tmp = UnivariateDivision.pseudoDivideAndRemainder(x, y, true);
            assert tmp != null && tmp[0] != null && tmp[1] != null;
            r = tmp[1];
            if (r.isZero())
                break;
            if (primitivePRS)
                r = r.primitivePart();
            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders<>(prs);
    }

    /**
     * Euclidean algorithm for polynomials over Euclidean rings that produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T>
    EuclidSubresultantRemainders(final T a, final T b) {
        if (a instanceof UnivariatePolynomialZ64)
            return (PolynomialRemainders<T>) EuclidSubresultantRemainders((UnivariatePolynomialZ64) a, (UnivariatePolynomialZ64) b);
        else if (a instanceof UnivariatePolynomial)
            return (PolynomialRemainders<T>) EuclidSubresultantRemainders((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        else
            throw new RuntimeException("Not a Z[x] polynomials: " + a.getClass());
    }

    /**
     * Euclidean algorithm for polynomials over Euclidean rings that produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    public static PolynomialRemainders<UnivariatePolynomialZ64>
    EuclidSubresultantRemainders(final UnivariatePolynomialZ64 a, final UnivariatePolynomialZ64 b) {
        if (b.degree > a.degree)
            return EuclidSubresultantRemainders(b, a);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());


        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        UnivariatePolynomialZ64 aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<UnivariatePolynomialZ64> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            UnivariatePolynomialZ64 curr = prs.get(i);
            UnivariatePolynomialZ64 next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? 1 : -1;
                cPsi = -1;
            } else {
                cPsi = MachineArithmetic.safePow(-curr.lc(), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = MachineArithmetic.safeMultiply(cPsi, MachineArithmetic.safePow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    long tmp = MachineArithmetic.safePow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    assert cPsi % tmp == 0;
                    cPsi /= tmp;
                }
                cBeta = MachineArithmetic.safeMultiply(-curr.lc(), MachineArithmetic.safePow(cPsi, delta));
            }

            UnivariatePolynomialZ64 q = UnivariateDivision.pseudoDivideAndRemainder(curr, next, true)[1];
            if (q.isZero())
                break;

            q = q.divideOrNull(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders<UnivariatePolynomialZ64> res = new PolynomialRemainders<>(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials over Euclidean rings that produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <E> PolynomialRemainders<UnivariatePolynomial<E>>
    EuclidSubresultantRemainders(final UnivariatePolynomial<E> a, final UnivariatePolynomial<E> b) {
        if (b.degree > a.degree)
            return EuclidSubresultantRemainders(b, a);

        if (a.isZero() || b.isZero())
            return new PolynomialRemainders<>(a.clone(), b.clone());


        Ring<E> ring = a.ring;
        E aContent = a.content(), bContent = b.content();
        E contentGCD = ring.gcd(aContent, bContent);
        UnivariatePolynomial<E> aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<UnivariatePolynomial<E>> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        ArrayList<E> beta = new ArrayList<>(), psi = new ArrayList<>();
        TIntArrayList deltas = new TIntArrayList();

        E cBeta, cPsi;
        for (int i = 0; ; i++) {
            UnivariatePolynomial<E> curr = prs.get(i);
            UnivariatePolynomial<E> next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? ring.getOne() : ring.getNegativeOne();
                cPsi = ring.getNegativeOne();
            } else {
                cPsi = ring.pow(ring.negate(curr.lc()), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = ring.multiply(cPsi, ring.pow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    E tmp = ring.pow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    //assert cPsi.remainder(tmp).isZero();
                    cPsi = ring.divideExact(cPsi, tmp);
                }
                cBeta = ring.multiply(ring.negate(curr.lc()), ring.pow(cPsi, delta));
            }

            UnivariatePolynomial<E> q = UnivariateDivision.pseudoDivideAndRemainder(curr, next, true)[1];
            if (q.isZero())
                break;

            q = q.divideOrNull(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders<UnivariatePolynomial<E>> res = new PolynomialRemainders<>(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }

    /**
     * Polynomial remainder sequence
     */
    public static final class PolynomialRemainders<T extends IUnivariatePolynomial<T>>
            implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        /** actual data */
        public final ArrayList<T> remainders;

        @SuppressWarnings("unchecked")
        public PolynomialRemainders(T... remainders) {
            this(new ArrayList<>(Arrays.asList(remainders)));
        }

        public PolynomialRemainders(ArrayList<T> remainders) {
            this.remainders = remainders;
        }

        public T gcd() {
            if (remainders.size() == 2 && remainders.get(1).isZero())
                return remainders.get(0);
            return remainders.get(remainders.size() - 1);
        }

        public int size() { return remainders.size(); }
    }

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static UnivariatePolynomialZ64 ModularGCD(UnivariatePolynomialZ64 a, UnivariatePolynomialZ64 b) {
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return UnivariatePolynomialZ64.create(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);
    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static UnivariatePolynomialZ64 ModularGCD0(UnivariatePolynomialZ64 a, UnivariatePolynomialZ64 b) {
        assert a.degree >= b.degree;

        long lcGCD = MachineArithmetic.gcd(a.lc(), b.lc());
        double bound = Math.max(a.mignotteBound(), b.mignotteBound()) * lcGCD;

        UnivariatePolynomialZ64 previousBase = null;
        UnivariatePolynomialZp64 base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            UnivariatePolynomialZp64 aMod = a.modulus(prime), bMod = b.modulus(prime);
            UnivariatePolynomialZp64 modularGCD = EuclidGCD(aMod, bMod);
            //clone if necessary
            if (modularGCD == aMod || modularGCD == bMod)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return UnivariatePolynomialZ64.one();

            // save the base for the first time or when a new modular image is better
            if (base == null || base.degree > modularGCD.degree) {
                //make base monic and multiply lcGCD
                modularGCD.monic(lcGCD);
                base = modularGCD;
                basePrime = prime;
                continue;
            }

            //skip unlucky prime
            if (base.degree < modularGCD.degree)
                continue;

            //lifting
            long newBasePrime = MachineArithmetic.safeMultiply(basePrime, prime);
            long monicFactor =
                    modularGCD.multiply(
                            MachineArithmetic.modInverse(modularGCD.lc(), prime),
                            modularGCD.ring.modulus(lcGCD));
            ChineseRemainders.ChineseRemaindersMagic magic = ChineseRemainders.createMagic(basePrime, prime);
            for (int i = 0; i <= base.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                //long oth = mod(safeMultiply(mod(safeMultiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);

                long oth = modularGCD.multiply(modularGCD.data[i], monicFactor);
                base.data[i] = ChineseRemainders(magic, base.data[i], oth);
            }
            base = base.setModulusUnsafe(newBasePrime);
            basePrime = newBasePrime;

            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            UnivariatePolynomialZ64 candidate = base.asPolyZSymmetric().primitivePart();
            if ((double) basePrime >= 2 * bound || (previousBase != null && candidate.equals(previousBase))) {
                previousBase = candidate;
                //first check b since b is less degree
                UnivariatePolynomialZ64[] div;
                div = UnivariateDivision.divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = UnivariateDivision.divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
            previousBase = candidate;
        }
    }

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("ConstantConditions")
    public static UnivariatePolynomial<BigInteger> ModularGCD(UnivariatePolynomial<BigInteger> a,
                                                              UnivariatePolynomial<BigInteger> b) {
        if (!a.ring.equals(Rings.Z))
            throw new IllegalArgumentException("Only polynomials over integers ring are allowed; " + a.ring);
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        BigInteger aContent = a.content(), bContent = b.content();
        BigInteger contentGCD = BigIntegerUtil.gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return a.createConstant(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);
    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static UnivariatePolynomial<BigInteger> ModularGCD0(UnivariatePolynomial<BigInteger> a,
                                                                UnivariatePolynomial<BigInteger> b) {
        assert a.degree >= b.degree;

        BigInteger lcGCD = BigIntegerUtil.gcd(a.lc(), b.lc());
        BigInteger bound2 = BigIntegerUtil.max(UnivariatePolynomial.mignotteBound(a), UnivariatePolynomial.mignotteBound(b)).multiply(lcGCD).shiftLeft(1);
        if (bound2.isLong()
                && UnivariatePolynomial.maxAbsCoefficient(a).isLong()
                && UnivariatePolynomial.maxAbsCoefficient(b).isLong())
            return ModularGCD(UnivariatePolynomial.asOverZ64(a), UnivariatePolynomial.asOverZ64(b)).toBigPoly();

        UnivariatePolynomialZ64 previousBase = null;
        UnivariatePolynomialZp64 base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(1031);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            BigInteger bPrime = BigInteger.valueOf(prime);
            if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                continue;

            IntegersZp bPrimeDomain = new IntegersZp(bPrime);
            UnivariatePolynomialZp64 aMod = UnivariatePolynomial.asOverZp64(a.setRing(bPrimeDomain)), bMod = UnivariatePolynomial.asOverZp64(b.setRing(bPrimeDomain));
            UnivariatePolynomialZp64 modularGCD = EuclidGCD(aMod, bMod);
            //clone if necessary
            if (modularGCD == aMod || modularGCD == bMod)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return a.createOne();

            // save the base for the first time or when a new modular image is better
            if (base == null || base.degree > modularGCD.degree) {
                //make base monic and multiply lcGCD
                long lLcGCD = lcGCD.mod(bPrime).longValueExact();
                modularGCD.monic(lLcGCD);
                base = modularGCD;
                basePrime = prime;
                continue;
            }

            //skip unlucky prime
            if (base.degree < modularGCD.degree)
                continue;

            if (MachineArithmetic.isOverflowMultiply(basePrime, prime) || basePrime * prime > MachineArithmetic.MAX_SUPPORTED_MODULUS)
                break;

            //lifting
            long newBasePrime = basePrime * prime;
            long monicFactor = modularGCD.multiply(
                    MachineArithmetic.modInverse(modularGCD.lc(), prime),
                    lcGCD.mod(bPrime).longValueExact());
            ChineseRemainders.ChineseRemaindersMagic magic = ChineseRemainders.createMagic(basePrime, prime);
            for (int i = 0; i <= base.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                //long oth = mod(safeMultiply(mod(safeMultiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);

                long oth = modularGCD.multiply(modularGCD.data[i], monicFactor);
                base.data[i] = ChineseRemainders(magic, base.data[i], oth);
            }
            base = base.setModulusUnsafe(newBasePrime);
            basePrime = newBasePrime;

            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            UnivariatePolynomialZ64 lCandidate = base.asPolyZSymmetric().primitivePart();
            if (BigInteger.valueOf(basePrime).compareTo(bound2) >= 0 || (previousBase != null && lCandidate.equals(previousBase))) {
                previousBase = lCandidate;
                UnivariatePolynomial<BigInteger> candidate = lCandidate.toBigPoly();
                //first check b since b is less degree
                UnivariatePolynomial<BigInteger>[] div;
                div = UnivariateDivision.divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = UnivariateDivision.divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
            previousBase = lCandidate;
        }

        //continue lifting with multi-precision integers
        UnivariatePolynomial<BigInteger> bPreviousBase = null, bBase = base.toBigPoly();
        BigInteger bBasePrime = BigInteger.valueOf(basePrime);

        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            BigInteger bPrime = BigInteger.valueOf(prime);
            if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                continue;

            IntegersZp bPrimeDomain = new IntegersZp(bPrime);
            UnivariatePolynomialZp64 aMod = UnivariatePolynomial.asOverZp64(a.setRing(bPrimeDomain)), bMod = UnivariatePolynomial.asOverZp64(b.setRing(bPrimeDomain));
            UnivariatePolynomialZp64 modularGCD = EuclidGCD(aMod, bMod);
            //clone if necessary
            if (modularGCD == aMod || modularGCD == bMod)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return a.createOne();

            //save the base
            if (bBase == null || bBase.degree > modularGCD.degree) {
                //make base monic and multiply lcGCD
                long lLcGCD = lcGCD.mod(bPrime).longValueExact();
                modularGCD.monic(lLcGCD);
                bBase = modularGCD.toBigPoly();
                bBasePrime = bPrime;
                continue;
            }

            //skip unlucky prime
            if (bBase.degree < modularGCD.degree)
                continue;

            //lifting
            BigInteger newBasePrime = bBasePrime.multiply(bPrime);
            long monicFactor = modularGCD.multiply(
                    MachineArithmetic.modInverse(modularGCD.lc(), prime),
                    lcGCD.mod(bPrime).longValueExact());
            for (int i = 0; i <= bBase.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                //long oth = mod(safeMultiply(mod(safeMultiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);

                long oth = modularGCD.multiply(modularGCD.data[i], monicFactor);
                bBase.data[i] = ChineseRemainders(bBasePrime, bPrime, bBase.data[i], BigInteger.valueOf(oth));
            }
            bBase = bBase.setRingUnsafe(new IntegersZp(newBasePrime));
            bBasePrime = newBasePrime;

            UnivariatePolynomial<BigInteger> candidate = UnivariatePolynomial.asPolyZSymmetric(bBase).primitivePart();
            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            if (bBasePrime.compareTo(bound2) >= 0 || (bPreviousBase != null && candidate.equals(bPreviousBase))) {
                bPreviousBase = candidate;
                //first check b since b is less degree
                UnivariatePolynomial<BigInteger>[] div;
                div = UnivariateDivision.divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = UnivariateDivision.divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
            bPreviousBase = candidate;
        }
    }
}
