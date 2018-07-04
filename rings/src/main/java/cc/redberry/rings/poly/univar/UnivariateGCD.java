package cc.redberry.rings.poly.univar;


import cc.redberry.rings.*;
import cc.redberry.rings.ChineseRemainders.ChineseRemaindersMagic;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.Util.Tuple2;
import cc.redberry.rings.poly.multivar.AMonomial;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateGCD;
import cc.redberry.rings.primes.PrimesIterator;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.ArraysUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static cc.redberry.rings.ChineseRemainders.ChineseRemainders;
import static cc.redberry.rings.ChineseRemainders.createMagic;
import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.Rings.Z;
import static cc.redberry.rings.poly.Util.toCommonDenominator;
import static cc.redberry.rings.poly.univar.Conversions64bit.asOverZp64;
import static cc.redberry.rings.poly.univar.Conversions64bit.canConvertToZp64;
import static cc.redberry.rings.poly.univar.UnivariatePolynomial.asOverZp64;

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
        if (a.isOverField()) {
            if (Util.isOverRationals(a))
                return (T) PolynomialGCDInQ((UnivariatePolynomial) a, (UnivariatePolynomial) b);
            return HalfGCD(a, b);
        } else if (a instanceof UnivariatePolynomialZ64)
            return (T) ModularGCD((UnivariatePolynomialZ64) a, (UnivariatePolynomialZ64) b);
        else if (a instanceof UnivariatePolynomial) {
            Ring ring = ((UnivariatePolynomial) a).ring;
            if (ring.equals(Z))
                return (T) ModularGCD((UnivariatePolynomial) a, (UnivariatePolynomial) b);

            if (ring instanceof AlgebraicNumberField) {
                Ring cfRing = ((UnivariatePolynomial) ((AlgebraicNumberField) ring).getMinimalPoly()).ring;
                if (cfRing.equals(Z))
                    return (T) PolynomialGCDInRingOfIntegersOfNumberField((UnivariatePolynomial) a, (UnivariatePolynomial) b);
                else if (cfRing.equals(Q))
                    return (T) PolynomialGCDInNumberField((UnivariatePolynomial) a, (UnivariatePolynomial) b);
            }

            T r = tryNested(a, b);
            if (r != null)
                return r;
            return (T) UnivariateResultants.SubresultantPRS((UnivariatePolynomial) a, (UnivariatePolynomial) b).gcd();
        } else
            throw new RuntimeException();
    }

    @SuppressWarnings("unchecked")
    private static <T extends IUnivariatePolynomial<T>> T tryNested(T a, T b) {
        if (a instanceof UnivariatePolynomial && ((UnivariatePolynomial) a).ring instanceof MultivariateRing)
            return (T) PolynomialGCDOverMultivariate((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        return null;
    }

    private static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    UnivariatePolynomial<Poly> PolynomialGCDOverMultivariate(UnivariatePolynomial<Poly> a,
                                                             UnivariatePolynomial<Poly> b) {
        return MultivariateGCD.PolynomialGCD(
                AMultivariatePolynomial.asMultivariate(a, 0, true),
                AMultivariatePolynomial.asMultivariate(b, 0, true))
                .asUnivariateEliminate(0);
    }

    private static <E> UnivariatePolynomial<Rational<E>> PolynomialGCDInQ(
            UnivariatePolynomial<Rational<E>> a,
            UnivariatePolynomial<Rational<E>> b) {

        Tuple2<UnivariatePolynomial<E>, E> aRat = toCommonDenominator(a);
        Tuple2<UnivariatePolynomial<E>, E> bRat = toCommonDenominator(b);

        return Util.asOverRationals(a.ring, PolynomialGCD(aRat._1, bRat._1)).monic();
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
        if (Util.isOverQ(a))
            return (T[]) ModularExtendedResultantGCDInQ((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        if (a.isOverZ())
            return (T[]) ModularExtendedResultantGCDInZ((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        if (a.isOverField())
            return ExtendedHalfGCD(a, b);
        else
            throw new IllegalArgumentException("Polynomial over field is expected");
    }

    /**
     * Returns array of {@code [gcd(a,b), s]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @param a the first poly for which the Bezout coefficient is computed
     * @param b the second poly
     * @return array of {@code [gcd(a,b), s]} such that {@code s * a + t * b = gcd(a, b)}
     */
    public static <T extends IUnivariatePolynomial<T>> T[] PolynomialFirstBezoutCoefficient(T a, T b) {
        if (a.isOverFiniteField() && Math.min(a.degree(), b.degree()) < 384)
            // this is somewhat faster than computing full xGCD
            return EuclidFirstBezoutCoefficient(a, b);
        else
            return Arrays.copyOf(PolynomialExtendedGCD(a, b), 2);
    }

    @SuppressWarnings("unchecked")
    static UnivariatePolynomial<BigInteger>[] PolynomialExtendedGCDInZbyQ(
            UnivariatePolynomial<BigInteger> a, UnivariatePolynomial<BigInteger> b) {
        UnivariatePolynomial<Rational<BigInteger>>[] xgcd = PolynomialExtendedGCD(
                a.mapCoefficients(Q, Q::mkNumerator),
                b.mapCoefficients(Q, Q::mkNumerator));

        Tuple2<UnivariatePolynomial<BigInteger>, BigInteger>
                gcd = toCommonDenominator(xgcd[0]),
                s = toCommonDenominator(xgcd[1]),
                t = toCommonDenominator(xgcd[2]);
        BigInteger lcm = Z.lcm(gcd._2, s._2, t._2);
        return new UnivariatePolynomial[]{
                gcd._1.multiply(lcm.divideExact(gcd._2)),
                s._1.multiply(lcm.divideExact(s._2)),
                t._1.multiply(lcm.divideExact(t._2))
        };
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

    static <T extends IUnivariatePolynomial<T>> T[] normalizeExtendedGCD(T[] xgcd) {
        if (!xgcd[0].isOverField())
            return xgcd;

        if (xgcd[0].isZero())
            return xgcd;
        xgcd[1].divideByLC(xgcd[0]);
        if (xgcd.length > 2)
            xgcd[2].divideByLC(xgcd[0]);
        xgcd[0].monic();
        return xgcd;
    }

    /**
     * Returns array of {@code [gcd(a,b), s]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @param a the first poly for which the Bezout coefficient is computed
     * @param b the second poly
     * @return array of {@code [gcd(a,b), s]} such that {@code s * a + t * b = gcd(a, b)}
     */
    public static <T extends IUnivariatePolynomial<T>> T[] EuclidFirstBezoutCoefficient(final T a, final T b) {
        a.assertSameCoefficientRingWith(b);
        if (canConvertToZp64(a))
            return Conversions64bit.convert(a, EuclidFirstBezoutCoefficient(asOverZp64(a), asOverZp64(b)));

        T s = a.createZero(), old_s = a.createOne();
        T r = b, old_r = a;

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
        }
        T[] result = a.createArray(2);
        result[0] = old_r;
        result[1] = old_s;
        return normalizeExtendedGCD(result);
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
            UnivariatePolynomialZp64 modularGCD = HalfGCD(aMod, bMod);
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
            ChineseRemainders.ChineseRemaindersMagicZp64 magic = ChineseRemainders.createMagic(basePrime, prime);
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
        if (!a.ring.equals(Z))
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
                && a.maxAbsCoefficient().isLong()
                && b.maxAbsCoefficient().isLong())
            try {
                // long overflow may occur here in very very rare cases
                return ModularGCD(UnivariatePolynomial.asOverZ64(a), UnivariatePolynomial.asOverZ64(b)).toBigPoly();
            } catch (ArithmeticException e) {}

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
            UnivariatePolynomialZp64 aMod = asOverZp64(a.setRing(bPrimeDomain)), bMod = asOverZp64(b.setRing(bPrimeDomain));
            UnivariatePolynomialZp64 modularGCD = HalfGCD(aMod, bMod);
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
            ChineseRemainders.ChineseRemaindersMagicZp64 magic = ChineseRemainders.createMagic(basePrime, prime);
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
            UnivariatePolynomialZp64 aMod = asOverZp64(a.setRing(bPrimeDomain)), bMod = asOverZp64(b.setRing(bPrimeDomain));
            UnivariatePolynomialZp64 modularGCD = HalfGCD(aMod, bMod);
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
            ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, bBasePrime, bPrime);
            for (int i = 0; i <= bBase.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                //long oth = mod(safeMultiply(mod(safeMultiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);

                long oth = modularGCD.multiply(modularGCD.data[i], monicFactor);
                bBase.data[i] = ChineseRemainders(Z, magic, bBase.data[i], BigInteger.valueOf(oth));
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

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings({"ConstantConditions", "unchecked"})
    public static UnivariatePolynomial<Rational<BigInteger>>[] ModularExtendedRationalGCD(
            UnivariatePolynomial<Rational<BigInteger>> a,
            UnivariatePolynomial<Rational<BigInteger>> b) {
        if (a == b || a.equals(b))
            return new UnivariatePolynomial[]{a.createOne(), a.createZero(), a.clone()};

        if (a.degree() < b.degree()) {
            UnivariatePolynomial<Rational<BigInteger>>[] r = ModularExtendedRationalGCD(b, a);
            ArraysUtil.swap(r, 1, 2);
            return r;
        }

        if (b.isZero()) {
            UnivariatePolynomial<Rational<BigInteger>>[] result = a.createArray(3);
            result[0] = a.clone();
            result[1] = a.createOne();
            result[2] = a.createZero();
            return normalizeExtendedGCD(result);
        }

        Tuple2<UnivariatePolynomial<BigInteger>, BigInteger>
                ac = toCommonDenominator(a),
                bc = toCommonDenominator(b);

        UnivariatePolynomial<BigInteger> az = ac._1, bz = bc._1;
        BigInteger
                aContent = az.content(),
                bContent = bz.content();

        UnivariatePolynomial<Rational<BigInteger>>[] xgcd = ModularExtendedRationalGCD0(
                az.clone().divideOrNull(aContent),
                bz.clone().divideOrNull(bContent));
        xgcd[1].multiply(new Rational<>(Z, ac._2, aContent));
        xgcd[2].multiply(new Rational<>(Z, bc._2, bContent));
        return xgcd;
    }

    /** modular extended GCD in Q[x] for primitive polynomials */
    @SuppressWarnings({"ConstantConditions", "unchecked"})
    static UnivariatePolynomial<Rational<BigInteger>>[] ModularExtendedRationalGCD0(
            UnivariatePolynomial<BigInteger> a,
            UnivariatePolynomial<BigInteger> b) {
        assert a.degree >= b.degree;
        BigInteger lcGCD = BigIntegerUtil.gcd(a.lc(), b.lc());

        UnivariatePolynomial<Rational<BigInteger>>
                aRat = a.mapCoefficients(Rings.Q, c -> new Rational<>(Z, c)),
                bRat = b.mapCoefficients(Rings.Q, c -> new Rational<>(Z, c));

        int degreeMax = Math.max(a.degree, b.degree);
        BigInteger bound2 = BigInteger.valueOf(degreeMax).increment().pow(degreeMax)
                .multiply(BigIntegerUtil.max(a.normMax(), b.normMax()).pow(a.degree + b.degree))
                .multiply(lcGCD)
                .shiftLeft(1);

        PrimesIterator primesLoop = new PrimesIterator(1031);//SmallPrimes.nextPrime(1 << 25));

        List<BigInteger> primes = new ArrayList<>();
        List<UnivariatePolynomial<BigInteger>>[] gst = new List[]{new ArrayList<>(), new ArrayList<>(), new ArrayList<>()};
        BigInteger primesMul = BigInteger.ONE;
        main:
        while (true) {
            while (primesMul.compareTo(bound2) < 0) {
                long prime = primesLoop.take();
                BigInteger bPrime = BigInteger.valueOf(prime);
                if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                    continue;

                IntegersZp bPrimeDomain = new IntegersZp(bPrime);
                UnivariatePolynomialZp64
                        aMod = asOverZp64(a.setRing(bPrimeDomain)),
                        bMod = asOverZp64(b.setRing(bPrimeDomain));

                UnivariatePolynomialZp64[] modularXGCD = ExtendedHalfGCD(aMod, bMod);
                assert modularXGCD[0].isMonic();

                //skip unlucky prime
                if (!gst[0].isEmpty() && modularXGCD[0].degree > gst[0].get(0).degree)
                    continue;

                if (!gst[0].isEmpty() && modularXGCD[0].degree < gst[0].get(0).degree) {
                    primes.clear();
                    primesMul = BigInteger.ONE;
                    Arrays.stream(gst).forEach(List::clear);
                }

                long lLcGCD = lcGCD.mod(bPrime).longValueExact();
                long lc = modularXGCD[0].lc();
                for (int i = 0; i < modularXGCD.length; i++)
                    gst[i].add(modularXGCD[i].multiply(lLcGCD).divide(lc).toBigPoly());
                primes.add(bPrime);
                primesMul = primesMul.multiply(bPrime);
            }


            // CRT
            UnivariatePolynomial<BigInteger>[] xgcdBase = new UnivariatePolynomial[3];
            BigInteger[] primesArray = primes.toArray(new BigInteger[primes.size()]);
            for (int i = 0; i < 3; ++i) {
                xgcdBase[i] = UnivariatePolynomial.zero(Z);
                int deg = gst[i].stream().mapToInt(UnivariatePolynomial::degree).max().getAsInt();
                xgcdBase[i].ensureCapacity(deg);
                for (int j = 0; j <= deg; ++j) {
                    final int jf = j;
                    BigInteger[] cfs = gst[i].stream().map(p -> p.get(jf)).toArray(BigInteger[]::new);
                    xgcdBase[i].data[j] = ChineseRemainders(primesArray, cfs);
                }
                xgcdBase[i].fixDegree();
            }

            while (true) {
                // do rational reconstruction
                UnivariatePolynomial<Rational<BigInteger>>[] xgcd = reconstructXGCD(aRat, bRat, xgcdBase, primesMul, bound2);
                if (xgcd != null)
                    return xgcd;

                // continue with CRT
                long prime = primesLoop.take();
                BigInteger bPrime = BigInteger.valueOf(prime);
                if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                    continue;

                IntegersZp bPrimeDomain = new IntegersZp(bPrime);
                UnivariatePolynomialZp64
                        aMod = asOverZp64(a.setRing(bPrimeDomain)),
                        bMod = asOverZp64(b.setRing(bPrimeDomain));

                UnivariatePolynomialZp64[] modularXGCD = ExtendedHalfGCD(aMod, bMod);
                assert modularXGCD[0].isMonic();

                //skip unlucky prime
                if (modularXGCD[0].degree > xgcdBase[0].degree)
                    continue;

                if (modularXGCD[0].degree < xgcdBase[0].degree) {
                    primes.clear();
                    Arrays.stream(gst).forEach(List::clear);

                    long lLcGCD = lcGCD.mod(bPrime).longValueExact();
                    long lc = modularXGCD[0].lc();
                    for (int i = 0; i < modularXGCD.length; i++)
                        gst[i].add(modularXGCD[i].multiply(lLcGCD).divide(lc).toBigPoly());
                    primes.add(bPrime);
                    primesMul = bPrime;

                    continue main; // <- extremely rare
                }

                long lLcGCD = lcGCD.mod(bPrime).longValueExact();
                long lc = modularXGCD[0].lc();
                for (UnivariatePolynomialZp64 m : modularXGCD)
                    m.multiply(lLcGCD).divide(lc);

                ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, primesMul, bPrime);
                for (int i = 0; i < 3; i++) {
                    modularXGCD[i].ensureCapacity(xgcdBase[i].degree);
                    assert modularXGCD[i].degree <= xgcdBase[i].degree;
                    for (int j = 0; j <= xgcdBase[i].degree; ++j)
                        xgcdBase[i].data[j] = ChineseRemainders(Z, magic, xgcdBase[i].data[j], BigInteger.valueOf(modularXGCD[i].data[j]));
                }
                primes.add(bPrime);
                primesMul = primesMul.multiply(bPrime);
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static UnivariatePolynomial<Rational<BigInteger>>[] reconstructXGCD(
            UnivariatePolynomial<Rational<BigInteger>> aRat, UnivariatePolynomial<Rational<BigInteger>> bRat,
            UnivariatePolynomial<BigInteger>[] xgcdBase, BigInteger prime, BigInteger bound2) {
        UnivariatePolynomial<Rational<BigInteger>>[] candidate = new UnivariatePolynomial[3];
        for (int i = 0; i < 3; i++) {
            candidate[i] = UnivariatePolynomial.zero(Rings.Q);
            candidate[i].ensureCapacity(xgcdBase[i].degree);
            for (int j = 0; j <= xgcdBase[i].degree; ++j) {
                BigInteger[] numDen = RationalReconstruction.reconstruct(xgcdBase[i].data[j], prime, bound2, bound2);
                if (numDen == null)
                    return null;
                candidate[i].data[j] = new Rational<>(Z, numDen[0], numDen[1]);
            }
            candidate[i].fixDegree();
        }

        BigInteger content = candidate[0].mapCoefficients(Z, Rational::numerator).content();
        Rational<BigInteger> corr = new Rational<>(Z, Z.getOne(), content);

        UnivariatePolynomial<Rational<BigInteger>>
                sCandidate = candidate[1].multiply(corr),
                tCandidate = candidate[2].multiply(corr),
                gCandidate = candidate[0].multiply(corr);

        //first check b since b is less degree
        UnivariatePolynomial<Rational<BigInteger>>[] bDiv;
        bDiv = UnivariateDivision.divideAndRemainder(bRat, gCandidate, true);
        if (!bDiv[1].isZero())
            return null;
        UnivariatePolynomial<Rational<BigInteger>>[] aDiv;
        aDiv = UnivariateDivision.divideAndRemainder(aRat, gCandidate, true);
        if (!aDiv[1].isZero())
            return null;

        if (!satisfiesXGCD(aDiv[0], sCandidate, bDiv[0], tCandidate))
            return null;
        return candidate;
    }

    @SuppressWarnings("unchecked")
    private static boolean satisfiesXGCD(
            UnivariatePolynomial<Rational<BigInteger>> a, UnivariatePolynomial<Rational<BigInteger>> s,
            UnivariatePolynomial<Rational<BigInteger>> b, UnivariatePolynomial<Rational<BigInteger>> t) {
        Rational<BigInteger>
                zero = Rational.zero(Z),
                one = Rational.one(Z);
        for (Rational<BigInteger> subs : new Rational[]{zero, one}) {
            Rational<BigInteger>
                    ea = a.evaluate(subs),
                    es = s.evaluate(subs),
                    eb = b.evaluate(subs),
                    et = t.evaluate(subs);

            if (!ea.multiply(es).add(eb.multiply(et)).isOne())
                return false;
        }
        return a.multiply(s).add(b.multiply(t)).isOne();
    }

    /**
     * Modular extended GCD algorithm for polynomials over Q with the use of resultants.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings({"ConstantConditions", "unchecked"})
    public static UnivariatePolynomial<Rational<BigInteger>>[]
    ModularExtendedResultantGCDInQ(UnivariatePolynomial<Rational<BigInteger>> a,
                                   UnivariatePolynomial<Rational<BigInteger>> b) {
        Tuple2<UnivariatePolynomial<BigInteger>, BigInteger>
                ra = toCommonDenominator(a),
                rb = toCommonDenominator(b);

        UnivariatePolynomial<BigInteger>[] xgcdZ = ModularExtendedResultantGCDInZ(ra._1, rb._1);
        BigInteger content = Z.gcd(xgcdZ[0].content(), ra._2, rb._2);
        xgcdZ[0].divideExact(content);

        UnivariatePolynomial<Rational<BigInteger>>[] xgcd =
                Arrays.stream(xgcdZ)
                        .map(p -> p.mapCoefficients(Q, Q::mkNumerator))
                        .toArray(UnivariatePolynomial[]::new);

        xgcd[1].multiply(Q.mkNumerator(ra._2.divideExact(content)));
        xgcd[2].multiply(Q.mkNumerator(rb._2.divideExact(content)));
        return xgcd;
    }

    /**
     * Modular extended GCD algorithm for polynomials over Z with the use of resultants.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings({"ConstantConditions", "unchecked"})
    public static UnivariatePolynomial<BigInteger>[]
    ModularExtendedResultantGCDInZ(UnivariatePolynomial<BigInteger> a,
                                   UnivariatePolynomial<BigInteger> b) {
        if (a == b || a.equals(b))
            return new UnivariatePolynomial[]{a.createOne(), a.createZero(), a.clone()};

        if (a.degree() < b.degree()) {
            UnivariatePolynomial<BigInteger>[] r = ModularExtendedResultantGCDInZ(b, a);
            ArraysUtil.swap(r, 1, 2);
            return r;
        }

        if (b.isZero()) {
            UnivariatePolynomial<BigInteger>[] result = a.createArray(3);
            result[0] = a.clone();
            result[1] = a.createOne();
            result[2] = a.createZero();
            return normalizeExtendedGCD(result);
        }

        BigInteger
                aContent = a.content(),
                bContent = b.content();

        a = a.clone().divideExact(aContent);
        b = b.clone().divideExact(bContent);

        UnivariatePolynomial<BigInteger> gcd = PolynomialGCD(a, b);
        a = UnivariateDivision.divideExact(a, gcd, false);
        b = UnivariateDivision.divideExact(b, gcd, false);

        UnivariatePolynomial<BigInteger>[] xgcd = ModularExtendedResultantGCD0(a, b);
        xgcd[0].multiply(gcd);

        UnivariatePolynomial<BigInteger> g = xgcd[0], s = xgcd[1], t = xgcd[2];

        BigInteger
                as = Z.gcd(aContent, s.content()),
                bt = Z.gcd(bContent, t.content());
        aContent = aContent.divideExact(as);
        bContent = bContent.divideExact(bt);

        s.divideExact(as);
        t.divideExact(bt);

        t.multiply(aContent);
        g.multiply(aContent);
        s.multiply(bContent);
        g.multiply(bContent);

        return xgcd;
    }

    /** modular extended GCD for primitive coprime polynomials */
    @SuppressWarnings({"ConstantConditions", "unchecked"})
    private static UnivariatePolynomial<BigInteger>[] ModularExtendedResultantGCD0(UnivariatePolynomial<BigInteger> a,
                                                                                   UnivariatePolynomial<BigInteger> b) {
        assert a.degree >= b.degree;

        BigInteger gcd = UnivariateResultants.ModularResultant(a, b);
        UnivariatePolynomial<BigInteger>[] previousBase = null, base = null;
        BigInteger basePrime = null;

        PrimesIterator primesLoop = new PrimesIterator(SmallPrimes.nextPrime(1 << 28));
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            BigInteger bPrime = BigInteger.valueOf(prime);
            if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                continue;

            IntegersZp ring = new IntegersZp(bPrime);
            UnivariatePolynomialZp64
                    aMod = asOverZp64(a.setRing(ring)),
                    bMod = asOverZp64(b.setRing(ring));

            UnivariatePolynomialZp64[] modularXGCD = PolynomialExtendedGCD(aMod, bMod);
            if (modularXGCD[0].degree != 0)
                continue;

            // xgcd over finite fields are always normalized
            assert modularXGCD[0].isOne();

            // resultant correction
            long correction = gcd.mod(bPrime).longValueExact();
            Arrays.stream(modularXGCD).forEach(p -> p.multiply(correction));

            //save the base
            if (base == null) {
                //make base monic and multiply lcGCD
                base = Arrays.stream(modularXGCD).map(UnivariatePolynomialZp64::toBigPoly).toArray(UnivariatePolynomial[]::new);
                basePrime = bPrime;
                continue;
            }

            //CRT lifting
            ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, basePrime, bPrime);
            BigInteger newBasePrime = basePrime.multiply(bPrime);
            for (int e = 0; e < 3; ++e) {
                base[e] = base[e].setRingUnsafe(new IntegersZp(newBasePrime));

                if (base[e].degree < modularXGCD[e].degree)
                    base[e].ensureCapacity(modularXGCD[e].degree);

                for (int i = 0; i <= base[e].degree; ++i)
                    base[e].data[i] = ChineseRemainders(Z, magic, base[e].get(i), BigInteger.valueOf(modularXGCD[e].get(i)));

                base[e].fixDegree();

            }
            basePrime = newBasePrime;

            // compute candidate
            UnivariatePolynomial<BigInteger>[] candidate = Arrays.stream(base)
                    .map(UnivariatePolynomial::asPolyZSymmetric)
                    .toArray(UnivariatePolynomial[]::new);
            BigInteger content = Z.gcd(candidate[0].content(), candidate[1].content(), candidate[2].content());
            Arrays.stream(candidate).forEach(p -> p.divideExact(content));
            // two trials didn't change the result, probably we are done
            if ((previousBase != null && Arrays.equals(candidate, previousBase))) {
                previousBase = candidate;
                if (!satisfiesXGCD(a, b, candidate))
                    continue;
                return candidate;
            }
            previousBase = candidate;
        }
    }

    @SuppressWarnings("unchecked")
    private static <E> boolean satisfiesXGCD(UnivariatePolynomial<E> a,
                                             UnivariatePolynomial<E> b,
                                             UnivariatePolynomial<E>[] xgcd) {
        Ring<E> ring = xgcd[0].ring;
        for (E subs : ring.createArray(ring.getZero(), ring.getOne())) {
            E
                    ea = a.evaluate(subs),
                    es = xgcd[1].evaluate(subs),
                    eb = b.evaluate(subs),
                    et = xgcd[2].evaluate(subs),
                    eg = xgcd[0].evaluate(subs);

            if (!ring.addMutable(ring.multiplyMutable(ea, es), ring.multiplyMutable(eb, et)).equals(eg))
                return false;
        }
        return a.clone().multiply(xgcd[1]).add(b.clone().multiply(xgcd[2])).equals(xgcd[0]);
    }

    ////////////////////////////////////// Modular GCD in algebraic number fields //////////////////////////////////////

    private static <E> UnivariatePolynomial<UnivariatePolynomial<E>>
    TrivialGCDInNumberField(UnivariatePolynomial<UnivariatePolynomial<E>> a, UnivariatePolynomial<UnivariatePolynomial<E>> b) {
        AlgebraicNumberField<UnivariatePolynomial<E>> ring
                = (AlgebraicNumberField<UnivariatePolynomial<E>>) a.ring;

        if (!a.stream().allMatch(ring::isSimpleNumber)
                || !b.stream().allMatch(ring::isSimpleNumber))
            return null;

        UnivariatePolynomial<E>
                ar = a.mapCoefficients(ring.getMinimalPoly().ring, UnivariatePolynomial::cc),
                br = b.mapCoefficients(ring.getMinimalPoly().ring, UnivariatePolynomial::cc);
        return PolynomialGCD(ar, br).mapCoefficients(ring, cf -> UnivariatePolynomial.constant(ring.getMinimalPoly().ring, cf));
    }

    /** Computes GCD via Langemyr & Mccallum modular algorithm over algebraic number field */
    @SuppressWarnings("ConstantConditions")
    public static UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
    PolynomialGCDInNumberField(UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a,
                               UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b) {
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> simpleGCD = TrivialGCDInNumberField(a, b);
        if (simpleGCD != null)
            return simpleGCD;

        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>>
                numberField = (AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>>) a.ring;
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly = numberField.getMinimalPoly();

        assert numberField.isField();

        a = a.clone();
        b = b.clone();

        // reduce problem to the case with integer monic minimal polynomial
        if (minimalPoly.stream().allMatch(Rational::isIntegral)) {
            // minimal poly is already monic & integer

            UnivariatePolynomial<BigInteger> minimalPolyZ = minimalPoly.mapCoefficients(Z, Rational::numerator);
            AlgebraicNumberField<UnivariatePolynomial<BigInteger>> numberFieldZ = new AlgebraicNumberField<>(minimalPolyZ);

            removeDenominators(a);
            removeDenominators(b);

            assert a.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral));
            assert b.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral));

            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcdZ =
                    gcdAssociateInNumberField(
                            a.mapCoefficients(numberFieldZ, cf -> cf.mapCoefficients(Z, Rational::numerator)),
                            b.mapCoefficients(numberFieldZ, cf -> cf.mapCoefficients(Z, Rational::numerator)));

            return gcdZ
                    .mapCoefficients(numberField, p -> p.mapCoefficients(Q, cf -> new Rational<>(Z, cf)))
                    .monic();
        } else {
            // replace s -> s / lc(minPoly)
            BigInteger minPolyLeadCoeff = toCommonDenominator(minimalPoly)._1.lc();
            Rational<BigInteger>
                    scale = new Rational<>(Z, Z.getOne(), minPolyLeadCoeff),
                    scaleReciprocal = scale.reciprocal();

            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>>
                    scaledNumberField = new AlgebraicNumberField<>(minimalPoly.scale(scale).monic());
            return PolynomialGCDInNumberField(
                    a.mapCoefficients(scaledNumberField, cf -> cf.scale(scale)),
                    b.mapCoefficients(scaledNumberField, cf -> cf.scale(scale)))
                    .mapCoefficients(numberField, cf -> cf.scale(scaleReciprocal));
        }
    }

    private static void pseudoMonicize(UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a) {
        UnivariatePolynomial<Rational<BigInteger>> inv = a.ring.reciprocal(a.lc());
        a.multiply(Util.toCommonDenominator(inv)._1.mapCoefficients(Q, Q::mkNumerator));
        assert a.lc().isConstant();
    }

    static BigInteger removeDenominators(UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a) {
        BigInteger denominator = Z.lcm(a.stream().map(p -> toCommonDenominator(p)._2).collect(Collectors.toList()));
        a.multiply(a.ring.valueOfBigInteger(denominator));
        return denominator;
    }

    /** Computes some GCD associate via Langemyr & Mccallum modular algorithm over algebraic integers */
    @SuppressWarnings("ConstantConditions")
    public static UnivariatePolynomial<UnivariatePolynomial<BigInteger>>
    PolynomialGCDInRingOfIntegersOfNumberField(UnivariatePolynomial<UnivariatePolynomial<BigInteger>> a,
                                               UnivariatePolynomial<UnivariatePolynomial<BigInteger>> b) {
        if (!a.lc().isConstant() || !b.lc().isConstant())
            throw new IllegalArgumentException("Univariate GCD in non-field extensions requires polynomials have integer leading coefficients.");

        UnivariatePolynomial<BigInteger>
                aContent = a.content(),
                bContent = b.content();
        assert aContent.isConstant();
        assert bContent.isConstant();

        UnivariatePolynomial<BigInteger> contentGCD = aContent.createConstant(aContent.cc().gcd(bContent.cc()));

        a = a.clone().divideExact(aContent);
        b = b.clone().divideExact(bContent);

        return gcdAssociateInNumberField0(a, b).multiply(contentGCD);
    }

    /** Computes some GCD associate via Langemyr & McCallum modular algorithm over algebraic integers */
    @SuppressWarnings("ConstantConditions")
    static UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcdAssociateInNumberField(
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> a,
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> b) {
        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> numberField
                = (AlgebraicNumberField<UnivariatePolynomial<BigInteger>>) a.ring;

        integerPrimitivePart(a);
        integerPrimitivePart(b);

        if (!a.lc().isConstant())
            a.multiply(numberField.cancellingMultiplier(a.lc()));

        if (!b.lc().isConstant())
            b.multiply(numberField.cancellingMultiplier(b.lc()));

        integerPrimitivePart(a);
        integerPrimitivePart(b);

        // if all coefficients are simple numbers (no algebraic elements)
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> simpleGCD = TrivialGCDInNumberField(a, b);
        if (simpleGCD != null)
            return simpleGCD;

        return gcdAssociateInNumberField0(a, b);
    }

    static BigInteger integerPrimitivePart(UnivariatePolynomial<UnivariatePolynomial<BigInteger>> p) {
        BigInteger gcd = Z.gcd(p.stream().flatMap(UnivariatePolynomial::stream).sorted().collect(Collectors.toList()));
        p.stream().forEach(cf -> cf.divideExact(gcd));
        return gcd;
    }

    /** Langemyr & McCallum modular algorithm for primitive polynomials with integer lead coefficients */
    @SuppressWarnings("ConstantConditions")
    static UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcdAssociateInNumberField0(
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> a,
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> b) {

        // conditions for Langemyr & McCallum algorithm
        assert a.lc().isConstant();
        assert b.lc().isConstant();

        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> numberField
                = (AlgebraicNumberField<UnivariatePolynomial<BigInteger>>) a.ring;
        UnivariateRing<UnivariatePolynomial<BigInteger>> auxRing
                = Rings.UnivariateRing(Z);

        UnivariatePolynomial<BigInteger> minimalPoly = numberField.getMinimalPoly();

        // Weinberger & Rothschild (1976) correction denominator
        BigInteger
                lcGCD = Z.gcd(a.lc().cc(), b.lc().cc()),
                disc = UnivariateResultants.Discriminant(minimalPoly),
                correctionFactor = disc.multiply(lcGCD);

        BigInteger crtPrime = null;
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcd = null, prevCandidate = null;
        PrimesIterator primes = new PrimesIterator(1 << 20);
        while (true) {
            long prime = primes.take();
            IntegersZp64 zpRing = new IntegersZp64(prime);
            UnivariatePolynomialZp64 minimalPolyMod = asOverZp64(minimalPoly, zpRing);
            if (minimalPolyMod.nNonZeroTerms() != minimalPoly.nNonZeroTerms())
                // bad prime
                continue;

            FiniteField<UnivariatePolynomialZp64> modRing = new FiniteField<>(minimalPolyMod);
            UnivariatePolynomial<UnivariatePolynomialZp64>
                    aMod = a.mapCoefficients(modRing, cf -> asOverZp64(cf, zpRing)),
                    bMod = b.mapCoefficients(modRing, cf -> asOverZp64(cf, zpRing));

            UnivariatePolynomial<UnivariatePolynomialZp64> gcdMod;
            try {
                gcdMod = PolynomialGCD(aMod, bMod);
            } catch (Throwable e) {
                // bad prime
                continue;
            }

            if (gcdMod.isConstant())
                return a.createOne();

            gcdMod.multiply(correctionFactor.mod(prime).longValue());

            BigInteger bPrime = BigInteger.valueOf(prime);
            if (crtPrime == null || gcdMod.degree < gcd.degree) {
                crtPrime = bPrime;
                gcd = gcdMod.mapCoefficients(auxRing, cf -> cf.toBigPoly().setRing(Z));
                continue;
            }

            if (gcdMod.degree > gcd.degree)
                // bad prime
                continue;

            ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, crtPrime, bPrime);
            boolean updated = false;
            for (int i = gcd.degree; i >= 0; --i) {
                boolean u = updateCRT(magic, gcd.data[i], gcdMod.data[i]);
                if (u) updated = true;
            }
            crtPrime = crtPrime.multiply(bPrime);

            // do trial division
            IntegersZp crtRing = new IntegersZp(crtPrime);
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> candidate =
                    gcd.mapCoefficients(numberField,
                            cf -> numberField.valueOf(UnivariatePolynomial.asPolyZSymmetric(cf.setRingUnsafe(crtRing))))
                            .primitivePart();

            if (prevCandidate == null) {
                prevCandidate = candidate;
                continue;
            }

            if (!updated || prevCandidate.equals(candidate)) {
                UnivariatePolynomial<UnivariatePolynomial<BigInteger>> rem;

                rem = UnivariateDivision.pseudoRemainderAdaptive(b, candidate, true);
                if (rem == null || !rem.isZero())
                    continue;

                rem = UnivariateDivision.pseudoRemainderAdaptive(a, candidate, true);
                if (rem == null || !rem.isZero())
                    continue;

                return candidate;
            }
            prevCandidate = candidate;
        }
    }

    /**
     * Apply CRT to a poly
     */
    public static boolean updateCRT(ChineseRemaindersMagic<BigInteger> magic,
                                    UnivariatePolynomial<BigInteger> accumulated,
                                    UnivariatePolynomialZp64 update) {
        boolean updated = false;
        accumulated.ensureCapacity(update.degree);
        for (int i = Math.max(accumulated.degree, update.degree); i >= 0; --i) {
            BigInteger oldCf = accumulated.get(i);
            BigInteger newCf = ChineseRemainders(Z, magic, oldCf, BigInteger.valueOf(update.get(i)));
            if (!oldCf.equals(newCf))
                updated = true;
            accumulated.data[i] = newCf;
        }
        accumulated.fixDegree();
        return updated;
    }
}
