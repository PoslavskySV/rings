package cc.r2.core.polynomial;


import cc.r2.core.number.primes.PrimesIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static cc.r2.core.polynomial.LongArithmetics.*;

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
        prs.add(a.clone().modulus(modulus)); prs.add(b.clone().modulus(modulus));

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
                cPsi = pow(-curr.lc(), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = multiply(cPsi, pow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    long tmp = pow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    assert cPsi % tmp == 0;
                    cPsi /= tmp;
                }
                cBeta = multiply(-curr.lc(), pow(cPsi, delta));
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

    public static MutableLongPoly ModularGCD(MutableLongPoly a, MutableLongPoly b) {
        if (a == b)
            return a.clone();
        if (a.degree < b.degree)
            return ModularGCD(b, a);

        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        return ModularGCD0(a.divide(aContent), b.divide(bContent)).multiply(contentGCD);

    }

    @SuppressWarnings("ConstantConditions")
    private static MutableLongPoly ModularGCD0(MutableLongPoly a, MutableLongPoly b) {
        if (a.degree < b.degree)
            return ModularGCD(b, a);

        long lcGCD = gcd(a.lc(), b.lc());
        double bound = Math.sqrt(a.degree + 1) * (1L << a.degree) * Math.max(a.norm(), b.norm()) * lcGCD;

        MutableLongPoly previousBase, base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            MutableLongPoly modularGCD = Euclid(a, b, prime).gcd();
            //clone if necessary
            if (modularGCD == a || modularGCD == b)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return MutableLongPoly.one();


            //save the base
            if (base == null) {
                //make base monic and multiply lcGCD
                modularGCD.monic(lcGCD, prime);
                base = modularGCD;
                basePrime = prime;
                continue;
            }

            //unlucky base => start over
            if (base.degree > modularGCD.degree) {
                base = null;
                basePrime = -1;
                continue;
            }

            //skip unlucky prime
            if (base.degree < modularGCD.degree)
                continue;

            //cache current base
            previousBase = base.clone();

            //lifting
            long newBasePrime = multiply(basePrime, prime);
            long monicFactor = modInverse(modularGCD.lc(), prime);
            long lcMod = mod(lcGCD, prime);
            for (int i = 0; i <= base.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                long oth = mod(multiply(mod(multiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);
                base.data[i] = ChineseRemainders(basePrime, prime, base.data[i], oth);
            }
            base.fixDegree();
            basePrime = newBasePrime;

            //we are lucky with Mignotte's bound
            if ((double) basePrime >= 2 * bound) {
                MutableLongPoly result = primitivePart((base.symModulus(basePrime)));
                assert pseudoDivideAndRemainderAdaptive(a, result)[1].isZero() && pseudoDivideAndRemainderAdaptive(b, result)[1].isZero();
                return result;
            }

            //two trials didn't change the result, probably we are done
            if (base.equals(previousBase)) {
                MutableLongPoly candidate = primitivePart(base.clone().symModulus(basePrime));
                //first check b since b is less degree
                if (pseudoDivideAndRemainderAdaptive(b, candidate)[1].isZero()
                        && pseudoDivideAndRemainderAdaptive(a, candidate)[1].isZero())
                    return candidate;
            }
        }
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
     * Divide and remainder
     *
     * @param dividend dividend
     * @param divider  divider
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider) {
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0) {
            MutableLongPoly div = dividend.clone().divideOrNull(divider.lc());
            if (div == null) return null;
            return new MutableLongPoly[]{div, MutableLongPoly.zero()};
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider);
        return divideAndRemainderGeneral0(dividend, divider, 1);
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
        if (dividend.degree < divider.degree)
            return null;
        long factor = pow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new MutableLongPoly[]{dividend.clone().multiply(factor / dividend.lc()), MutableLongPoly.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor);
        return divideAndRemainderGeneral0(dividend, divider, factor);
    }

    /** Plain school implementation */
    static MutableLongPoly[] divideAndRemainderGeneral0(final MutableLongPoly dividend,
                                                        final MutableLongPoly divider,
                                                        long dividendRaiseFactor) {
        assert dividend.degree >= divider.degree;

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
     * Pseudo divide and remainder
     *
     * @param dividend dividend
     * @param divider  divider
     * @return {quotient, remainder}
     */
    static MutableLongPoly[] pseudoDivideAndRemainderAdaptive(final MutableLongPoly dividend,
                                                              final MutableLongPoly divider) {
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0)
            return new MutableLongPoly[]{dividend.clone(), MutableLongPoly.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider);
    }

    /** general implementation */
    static MutableLongPoly[] pseudoDivideAndRemainderAdaptive0(
            final MutableLongPoly dividend,
            final MutableLongPoly divider) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly
                remainder = dividend.clone(),
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0) {
                    long gcd = gcd(remainder.lc(), divider.lc());
                    long factor = divider.lc() / gcd;
                    remainder.multiply(factor);
                    quotient.multiply(factor);
                }

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
     * @param dividend dividend
     * @param divider  divider
     * @param modulus  modulus
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       final long modulus) {
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0)
            return new MutableLongPoly[]{dividend.clone().multiply(modInverse(divider.lc(), modulus), modulus), MutableLongPoly.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, modulus);
        return divideAndRemainderModulus(dividend, divider, modulus);
    }

    /** Plain school implementation */
    static MutableLongPoly[] divideAndRemainderModulus(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       final long modulus) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly
                remainder = dividend.clone(),
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient.data[i] = divideMod(remainder.lc(), divider.lc(), modulus);
                remainder.subtract(divider, quotient.data[i], i, modulus);
            } else quotient.data[i] = 0;
        }

        return new MutableLongPoly[]{quotient.modulus(modulus), remainder.modulus(modulus)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDividerModulus(MutableLongPoly dividend, MutableLongPoly divider, long modulus) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = mod(-divider.cc(), modulus);
        long lcInverse = modInverse(divider.lc(), modulus);

        if (divider.lc() != 1)
            cc = mod(multiply(cc, lcInverse), modulus);

        long[] quotient = new long[dividend.degree];
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            if (i != dividend.degree)
                quotient[i] = mod(multiply(res, lcInverse), modulus);
            res = addMod(multiply(res, cc), dividend.data[i], modulus);
        }
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDivider(MutableLongPoly dividend, MutableLongPoly divider) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] pseudoDivideAndRemainderLinearDivider(MutableLongPoly dividend, MutableLongPoly divider) {
        return divideAndRemainderLinearDivider0(dividend, divider, pow(divider.lc(), dividend.degree));
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDivider0(MutableLongPoly dividend, MutableLongPoly divider, long raiseFactor) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc();
        long[] quotient = new long[dividend.degree];
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(raiseFactor, dividend.data[i]));
            if (i == 0) break;
            if (res % lc != 0) return null;
            res = res / lc;
        }
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] pseudoDivideAndRemainderLinearDividerAdaptive(MutableLongPoly dividend, MutableLongPoly divider) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc(), factor = 1;
        long[] quotient = new long[dividend.degree];
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(factor, dividend.data[i]));
            if (i == 0) break;
            if (res % lc != 0) {
                long gcd = gcd(res, lc), f = lc / gcd;
                factor = multiply(factor, f);
                res = multiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = multiply(quotient[j], f);
            }
            res = res / lc;
        }
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }

    /**
     * Representation for polynomial remainders sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders {
        public final ArrayList<MutableLongPoly> remainders;

        public PolynomialRemainders(ArrayList<MutableLongPoly> remainders) {
            this.remainders = remainders;
        }

        public MutableLongPoly gcd() {
            return remainders.get(remainders.size() - 1);
        }
    }
}
