package cc.r2.core.poly.univar;


import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.number.ChineseRemainders;
import cc.r2.core.number.primes.PrimesIterator;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Integers;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.MachineArithmetic;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static cc.r2.core.poly.MachineArithmetic.safeMultiply;
import static cc.r2.core.poly.MachineArithmetic.safePow;
import static cc.r2.core.poly.univar.DivisionWithRemainder.*;
import static cc.r2.core.poly.univar.UnivariatePolynomial.*;

/**
 * Polynomial GCD and sub-resultant sequence for univariate polynomials.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariateGCD {
    private UnivariateGCD() {}

    /**
     * Euclidean algorithm
     *
     * @param a poly
     * @param b poly
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> EuclidRemainders(final T a, final T b) {
        if (a.degree() < b.degree())
            return EuclidRemainders(b, a);

        ArrayList<T> prs = new ArrayList<>();
        prs.add(a.clone()); prs.add(b.clone());

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(prs);

        T x = a, y = b, r;
        while (true) {
            r = remainder(x, y, true);
            if (r == null)
                throw new IllegalArgumentException("Not divisible: (" + x + ") / (" + y + ")");

            if (r.isZero())
                break;
            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders<>(prs);
    }

    /**
     * Euclidean algorithm
     *
     * @param a poly
     * @param b poly
     * @return the GCD (monic if a and b are over field)
     */
    public static <T extends IUnivariatePolynomial<T>> T EuclidGCD(final T a, final T b) {
        if (a.degree() < b.degree())
            return EuclidGCD(b, a);

        if (a.isZero()) return normalizeGCD(b.clone());
        if (b.isZero()) return normalizeGCD(a.clone());

        T x = a, y = b, r;
        while (true) {
            r = remainder(x, y, true);
            if (r == null)
                throw new IllegalArgumentException("Not divisible: (" + x + ") / (" + y + ")");

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
     * Runs extended Euclidean algorithm to compute {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> T[] ExtendedEuclidGCD(final T a, final T b) {
        T s = a.createZero(), old_s = a.createOne();
        T t = a.createOne(), old_t = a.createZero();
        T r = b.clone(), old_r = a.clone();

        T q;
        T tmp;
        while (!r.isZero()) {
            q = quotient(old_r, r, true);
            if (q == null)
                throw new IllegalArgumentException("Not divisible: (" + old_r + ") / (" + r + ")");

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
     * Half-GCD algorithm. The algorithm automatically switches to Euclidean.
     *
     * @param a poly
     * @param b poly
     * @return the GCD (monic)
     */
    public static <T extends IUnivariatePolynomial<T>> T HalfGCD(T a, T b) {
        if (a.degree() < b.degree())
            return HalfGCD(b, a);

        if (a.degree() == b.degree())
            b = remainder(b, a, true);

        while (a.degree() > SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE && !b.isZero()) {
            T[] col = reduceHalfGCD(a, b);
            a = col[0];
            b = col[1];

            if (!b.isZero()) {
                T remainder = remainder(a, b, true);
                a = b;
                b = remainder;
            }
        }
        return UnivariateGCD.EuclidGCD(a, b);
    }

    /**
     * Runs extended Half-GCD algorithm to compute {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)} (gcd is monic)
     */
    public static <T extends IUnivariatePolynomial<T>> T[] ExtendedHalfGCD(T a, T b) {
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
            T[] qd = divideAndRemainder(a, b, true);
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
            T[] qd = divideAndRemainder(tmpA, tmpB, true);
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

        T[] qd = divideAndRemainder(a1, b1, false);
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

        T[] qd = divideAndRemainder(a, b, true);
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

        T remainder = remainder(a, b, true);
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
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> PseudoEuclid(final T a,
                                                                                            final T b,
                                                                                            boolean primitivePRS) {
        if (a instanceof lUnivariatePolynomialZ)
            return (PolynomialRemainders<T>) PseudoEuclid((lUnivariatePolynomialZ) a, (lUnivariatePolynomialZ) b, primitivePRS);
        else if (a instanceof UnivariatePolynomial)
            return (PolynomialRemainders<T>) PseudoEuclid((UnivariatePolynomial) a, (UnivariatePolynomial) b, primitivePRS);
        else
            throw new RuntimeException("Not a Z[x] polynomials: " + a.getClass());
    }

    /**
     * Euclidean algorithm for polynomials over Z that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to build primitive polynomial remainders or not
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static PolynomialRemainders<lUnivariatePolynomialZ> PseudoEuclid(final lUnivariatePolynomialZ a,
                                                                            final lUnivariatePolynomialZ b,
                                                                            boolean primitivePRS) {
        if (a.degree < b.degree)
            return PseudoEuclid(b, a, primitivePRS);


        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());

        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        lUnivariatePolynomialZ aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);
        PolynomialRemainders<lUnivariatePolynomialZ> res = PseudoEuclid0(aPP, bPP, primitivePRS);
        res.gcd().primitivePartSameSign().multiply(contentGCD);
        return res;
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
    public static <E> PolynomialRemainders<UnivariatePolynomial<E>> PseudoEuclid(final UnivariatePolynomial<E> a,
                                                                                 final UnivariatePolynomial<E> b,
                                                                                 boolean primitivePRS) {
        if (a.degree < b.degree)
            return PseudoEuclid(b, a, primitivePRS);

        if (a.isZero() || b.isZero())
            return new PolynomialRemainders<>(a.clone(), b.clone());

        E aContent = a.content(), bContent = b.content();
        E contentGCD = a.domain.gcd(aContent, bContent);
        UnivariatePolynomial<E> aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);
        PolynomialRemainders<UnivariatePolynomial<E>> res = PseudoEuclid0(aPP, bPP, primitivePRS);
        res.gcd().primitivePartSameSign().multiply(contentGCD);
        return res;
    }

    @SuppressWarnings("unchecked")
    private static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> PseudoEuclid0(final T aPP,
                                                                                              final T bPP,
                                                                                              boolean primitivePRS) {
        ArrayList<T> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        T x = aPP, y = bPP, r;
        while (true) {
            T[] tmp = DivisionWithRemainder.pseudoDivideAndRemainder(x, y, true);
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
     * Euclidean algorithm for polynomials over Z that builds subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> SubresultantEuclid(final T a,
                                                                                                  final T b) {
        if (a instanceof lUnivariatePolynomialZ)
            return (PolynomialRemainders<T>) SubresultantEuclid((lUnivariatePolynomialZ) a, (lUnivariatePolynomialZ) b);
        else if (a instanceof UnivariatePolynomial)
            return (PolynomialRemainders<T>) SubresultantEuclid((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        else
            throw new RuntimeException("Not a Z[x] polynomials: " + a.getClass());
    }

    /**
     * Euclidean algorithm for polynomials over Z that builds subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    public static PolynomialRemainders<lUnivariatePolynomialZ> SubresultantEuclid(final lUnivariatePolynomialZ a,
                                                                                  final lUnivariatePolynomialZ b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());


        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        lUnivariatePolynomialZ aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<lUnivariatePolynomialZ> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            lUnivariatePolynomialZ curr = prs.get(i);
            lUnivariatePolynomialZ next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? 1 : -1;
                cPsi = -1;
            } else {
                cPsi = safePow(-curr.lc(), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = safeMultiply(cPsi, safePow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    long tmp = safePow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    assert cPsi % tmp == 0;
                    cPsi /= tmp;
                }
                cBeta = safeMultiply(-curr.lc(), safePow(cPsi, delta));
            }

            lUnivariatePolynomialZ q = DivisionWithRemainder.pseudoDivideAndRemainder(curr, next, true)[1];
            if (q.isZero())
                break;

            q = q.divideOrNull(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders<lUnivariatePolynomialZ> res = new PolynomialRemainders<>(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials over Z that builds subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    @SuppressWarnings("unchecked")
    public static <E> PolynomialRemainders<UnivariatePolynomial<E>> SubresultantEuclid(final UnivariatePolynomial<E> a,
                                                                                       final UnivariatePolynomial<E> b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        if (a.isZero() || b.isZero())
            return new PolynomialRemainders<>(a.clone(), b.clone());


        Domain<E> domain = a.domain;
        E aContent = a.content(), bContent = b.content();
        E contentGCD = domain.gcd(aContent, bContent);
        UnivariatePolynomial<E> aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<UnivariatePolynomial<E>> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        ArrayList<E> beta = new ArrayList<>(), psi = new ArrayList<>();
        TIntArrayList deltas = new TIntArrayList();

        E cBeta, cPsi;
        for (int i = 0; ; i++) {
            UnivariatePolynomial<E> curr = prs.get(i);
            UnivariatePolynomial<E> next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? domain.getOne() : domain.getNegativeOne();
                cPsi = domain.getNegativeOne();
            } else {
                cPsi = domain.pow(domain.negate(curr.lc()), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = domain.multiply(cPsi, domain.pow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    E tmp = domain.pow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    //assert cPsi.remainder(tmp).isZero();
                    cPsi = domain.divideExact(cPsi, tmp);
                }
                cBeta = domain.multiply(domain.negate(curr.lc()), domain.pow(cPsi, delta));
            }

            UnivariatePolynomial<E> q = pseudoDivideAndRemainder(curr, next, true)[1];
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
     * Polynomial remainder sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders<T extends IUnivariatePolynomial<T>> {
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
    public static lUnivariatePolynomialZ ModularGCD(lUnivariatePolynomialZ a, lUnivariatePolynomialZ b) {
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        long aContent = a.content(), bContent = b.content();
        long contentGCD = MachineArithmetic.gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return lUnivariatePolynomialZ.create(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);
    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static lUnivariatePolynomialZ ModularGCD0(lUnivariatePolynomialZ a, lUnivariatePolynomialZ b) {
        assert a.degree >= b.degree;

        long lcGCD = MachineArithmetic.gcd(a.lc(), b.lc());
        double bound = Math.max(a.mignotteBound(), b.mignotteBound()) * lcGCD;

        lUnivariatePolynomialZ previousBase = null;
        lUnivariatePolynomialZp base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            lUnivariatePolynomialZp aMod = a.modulus(prime), bMod = b.modulus(prime);
            lUnivariatePolynomialZp modularGCD = EuclidGCD(aMod, bMod);
            //clone if necessary
            if (modularGCD == aMod || modularGCD == bMod)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return lUnivariatePolynomialZ.one();

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
            long newBasePrime = safeMultiply(basePrime, prime);
            long monicFactor =
                    modularGCD.multiply(
                            MachineArithmetic.modInverse(modularGCD.lc(), prime),
                            modularGCD.domain.modulus(lcGCD));
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
            lUnivariatePolynomialZ candidate = base.asPolyZSymmetric().primitivePart();
            if ((double) basePrime >= 2 * bound || (previousBase != null && candidate.equals(previousBase))) {
                previousBase = candidate;
                //first check b since b is less degree
                lUnivariatePolynomialZ[] div;
                div = DivisionWithRemainder.divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = DivisionWithRemainder.divideAndRemainder(a, candidate, true);
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
        if (a.domain != Integers.Integers)
            throw new IllegalArgumentException("Only polynomials over integers domain are allowed; " + a.domain);
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        BigInteger aContent = a.content(), bContent = b.content();
        BigInteger contentGCD = BigIntegerArithmetics.gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return a.createConstant(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);
    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static UnivariatePolynomial<BigInteger> ModularGCD0(UnivariatePolynomial<BigInteger> a,
                                                                UnivariatePolynomial<BigInteger> b) {
        assert a.degree >= b.degree;

        BigInteger lcGCD = BigIntegerArithmetics.gcd(a.lc(), b.lc());
        BigInteger bound2 = BigIntegerArithmetics.max(mignotteBound(a), mignotteBound(b)).multiply(lcGCD).shiftLeft(1);
        if (bound2.isLong()
                && maxAbsCoefficient(a).isLong()
                && maxAbsCoefficient(b).isLong())
            return ModularGCD(asLongPolyZ(a), asLongPolyZ(b)).toBigPoly();

        lUnivariatePolynomialZ previousBase = null;
        lUnivariatePolynomialZp base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(1031);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            BigInteger bPrime = BigInteger.valueOf(prime);
            if (a.lc().remainder(bPrime).isZero() || b.lc().remainder(bPrime).isZero())
                continue;

            IntegersModulo bPrimeDomain = new IntegersModulo(bPrime);
            lUnivariatePolynomialZp aMod = asLongPolyZp(a.setDomain(bPrimeDomain)), bMod = asLongPolyZp(b.setDomain(bPrimeDomain));
            lUnivariatePolynomialZp modularGCD = EuclidGCD(aMod, bMod);
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
            lUnivariatePolynomialZ lCandidate = base.asPolyZSymmetric().primitivePart();
            if (BigInteger.valueOf(basePrime).compareTo(bound2) >= 0 || (previousBase != null && lCandidate.equals(previousBase))) {
                previousBase = lCandidate;
                UnivariatePolynomial<BigInteger> candidate = lCandidate.toBigPoly();
                //first check b since b is less degree
                UnivariatePolynomial<BigInteger>[] div;
                div = divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = divideAndRemainder(a, candidate, true);
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

            IntegersModulo bPrimeDomain = new IntegersModulo(bPrime);
            lUnivariatePolynomialZp aMod = asLongPolyZp(a.setDomain(bPrimeDomain)), bMod = asLongPolyZp(b.setDomain(bPrimeDomain));
            lUnivariatePolynomialZp modularGCD = EuclidGCD(aMod, bMod);
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
            bBase = bBase.setDomainUnsafe(new IntegersModulo(newBasePrime));
            bBasePrime = newBasePrime;

            UnivariatePolynomial<BigInteger> candidate = asPolyZSymmetric(bBase).primitivePart();
            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            if (bBasePrime.compareTo(bound2) >= 0 || (bPreviousBase != null && candidate.equals(bPreviousBase))) {
                bPreviousBase = candidate;
                //first check b since b is less degree
                UnivariatePolynomial<BigInteger>[] div;
                div = divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
            bPreviousBase = candidate;
        }
    }

    /**
     * Returns GCD of two polynomials. Modular GCD algorithm is used for Z[x], plain Euclid is used for Zp[x] and
     * subresultant used for other domains.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> T PolynomialGCD(T a, T b) {
        if (a.isOverField())
            return HalfGCD(a, b);

        if (a instanceof lUnivariatePolynomialZ)
            return (T) ModularGCD((lUnivariatePolynomialZ) a, (lUnivariatePolynomialZ) b);
        else if (a instanceof UnivariatePolynomial) {
            Domain domain = ((UnivariatePolynomial) a).domain;
            if (domain == Integers.Integers)
                return (T) ModularGCD((UnivariatePolynomial) a, (UnivariatePolynomial) b);
            else
                return (T) SubresultantEuclid((UnivariatePolynomial) a, (UnivariatePolynomial) b).gcd();
        } else
            throw new RuntimeException(a.getClass().toString());
    }

    /**
     * Computes {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)} (gcd is monic)
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
}
