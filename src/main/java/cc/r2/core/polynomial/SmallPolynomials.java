package cc.r2.core.polynomial;


import cc.r2.core.number.ArithmeticUtils;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

import static cc.r2.core.number.ArithmeticUtils.*;

public final class SmallPolynomials {
    private SmallPolynomials() {
    }

    /**
     * Plain Euclidean algorithm which fails if intermediate polynomials are not divisible at some step
     *
     * @param a poly
     * @param b poly
     * @return a list of polynomial remainders where the last element is GCD
     * @throws IllegalArgumentException if at some step intermediate polynomials are not divisible
     */
    public static PolynomialRemainders Euclid(final MutableLongPoly a,
                                              final MutableLongPoly b) {
        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(a); prs.add(b);

        MutableLongPoly x = a, y = b, r;
        while (true) {
            MutableLongPoly[] tmp = divideAndRemainder(x, y);
            if (tmp == null)
                throw new IllegalArgumentException("Not divisible: (" + x + ") / (" + y + ")");

            r = tmp[1];
            if (r.isZero())
                break;

            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders(prs);
    }

    /**
     * Euclidean algorithm
     *
     * @param a poly
     * @param b poly
     * @return a list of polynomial remainders where the last element is GCD
     */
    public static PolynomialRemainders Euclid(final MutableLongPoly a,
                                              final MutableLongPoly b,
                                              final long modulus) {
        if (a.degree < b.degree)
            return Euclid(b, a, modulus);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(a); prs.add(b);

        MutableLongPoly x = a, y = b, r;
        while (true) {
            MutableLongPoly[] tmp = divideAndRemainder(x, y, modulus);
            assert tmp != null;
            r = tmp[1];
            if (r.isZero())
                break;
            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders(prs);
    }

    /**
     * Euclidean algorithm for polynomials that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to use primitive polynomial remainders
     * @return a list of polynomial remainders where the last element is GCD
     */
    public static PolynomialRemainders PolynomialEuclid(final MutableLongPoly a,
                                                        final MutableLongPoly b,
                                                        boolean primitivePRS) {
        if (a.degree < b.degree)
            return PolynomialEuclid(b, a, primitivePRS);

        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        MutableLongPoly aPP = a.divide(aContent), bPP = b.divide(bContent);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        MutableLongPoly x = aPP, y = bPP, r;
        while (true) {
            MutableLongPoly[] tmp = pseudoDivideAndRemainder(x, y);
            assert tmp != null;
            r = tmp[1];
            if (r.isZero())
                break;
            if (primitivePRS)
                r = primitivePart(r);
            prs.add(r);
            x = y;
            y = r;
        }
        PolynomialRemainders res = new PolynomialRemainders(prs);
        primitivePart(res.gcd()).multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials which produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence where the last element is GCD
     */
    public static PolynomialRemainders SubresultantEuclid(final MutableLongPoly a,
                                                          final MutableLongPoly b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        MutableLongPoly aPP = a.divide(aContent), bPP = b.divide(bContent);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            MutableLongPoly curr = prs.get(i);
            MutableLongPoly next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? 1 : -1;
                cPsi = -1;
            } else {
                assert safePow(-curr.lc(), deltas.get(i - 1));
                cPsi = pow(-curr.lc(), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    assert safePow(psi.get(i - 1), -deltas.get(i - 1) + 1);
                    assert safeMultiply(cPsi, pow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                    cPsi = cPsi * pow(psi.get(i - 1), -deltas.get(i - 1) + 1);
                } else {
                    assert safePow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    long tmp = pow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    assert cPsi % tmp == 0;
                    cPsi /= tmp;
                }
                assert safePow(cPsi, delta);
                assert safeMultiply(-curr.lc(), pow(cPsi, delta));
                cBeta = -curr.lc() * pow(cPsi, delta);
            }

            MutableLongPoly q = pseudoDivideAndRemainder(curr, next)[1];
            if (q.isZero())
                break;

            q = q.divide(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders res = new PolynomialRemainders(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }

    /**
     * Returns the content of the poly
     *
     * @param poly poly
     * @return polynomial content
     */
    public static long content(MutableLongPoly poly) {
        if (poly.degree == 0)
            return poly.data[0];
        return gcd(poly.data, 0, poly.degree + 1);
    }

    /**
     * Reduces poly to its primitive part
     *
     * @param poly polynomial
     * @return primitive part (poly will be modified)
     */
    public static MutableLongPoly primitivePart(MutableLongPoly poly) {
        long content = content(poly);
        if (content == 1)
            return poly;
        if (poly.lc() < 0)
            content = -content;
        for (int i = 0; i <= poly.degree; ++i) {
            assert poly.data[i] % content == 0;
            poly.data[i] = poly.data[i] / content;
        }
        return poly;
    }

    /**
     * Pseudo divide and remainder
     *
     * @param dividend dividend
     * @param divider  divider
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] pseudoDivideAndRemainder(
            MutableLongPoly dividend,
            MutableLongPoly divider) {
        assert safePow(divider.lc(), dividend.degree - divider.degree + 1) : "long overflow";
        return divideAndRemainder0(dividend, divider, pow(divider.lc(), dividend.degree - divider.degree + 1));
    }

    /**
     * Divide and remainder
     *
     * @param dividend dividend
     * @param divider  divider
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider) {
        return divideAndRemainder0(dividend, divider, 1);
    }

    /**
     * Divide and remainder
     *
     * @param dividend dividend
     * @param divider  divider
     * @return {quotient, remainder}
     */
    private static MutableLongPoly[] divideAndRemainder0(final MutableLongPoly dividend,
                                                         final MutableLongPoly divider,
                                                         long dividendRaiseFactor) {
        if (dividend.degree < divider.degree)
            return null;

        if (divider.degree == 0) {
            if (dividendRaiseFactor % divider.lc() == 0)
                return new MutableLongPoly[]{dividend.clone().multiply(dividendRaiseFactor / divider.lc()), MutableLongPoly.zero()};
            else {
                MutableLongPoly quot = dividend.clone().multiply(dividendRaiseFactor).divideOrNull(divider.lc());
                if (quot == null)
                    return null;
                return new MutableLongPoly[]{quot, MutableLongPoly.zero()};
            }
        }

        MutableLongPoly
                remainder = dividend.clone().multiply(dividendRaiseFactor),
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0)
                    return null;

                quotient.data[i] = remainder.lc() / divider.lc();
                remainder.subtract(divider, quotient.data[i], i);

            } else quotient.data[i] = 0;
        }

        quotient.fixDegree();
        return new MutableLongPoly[]{quotient, remainder};
    }

    /**
     * Divide and remainder
     *
     * @param dividend divident
     * @param divider  divider
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       final long modulus) {
        if (dividend.degree < divider.degree)
            return null;

        MutableLongPoly
                remainder = dividend.clone(),
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient.data[i] = ArithmeticUtils.divide(remainder.lc(), divider.lc(), modulus);
                remainder.subtract(divider, quotient.data[i], i, modulus);
            } else quotient.data[i] = 0;
        }

        return new MutableLongPoly[]{quotient.modulus(modulus), remainder.modulus(modulus)};
    }

    /**
     * Representation for polynomial remainders sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders {
        public final ArrayList<MutableLongPoly> remainders;

        public PolynomialRemainders(ArrayList<MutableLongPoly> remainders) {
            this.remainders = remainders;
//            if (remainders.size() == 2)
//                remainders.add(MutableLongPoly.one());
        }

        public MutableLongPoly gcd() {
            return remainders.get(remainders.size() - 1);
        }
    }
}
