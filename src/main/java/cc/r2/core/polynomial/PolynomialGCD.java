package cc.r2.core.polynomial;


import cc.r2.core.number.primes.PrimesIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static cc.r2.core.polynomial.DivideAndRemainder.divideAndRemainder;
import static cc.r2.core.polynomial.DivideAndRemainder.pseudoDivideAndRemainder;
import static cc.r2.core.polynomial.LongArithmetics.*;

public final class PolynomialGCD {
    private PolynomialGCD() {}

    /**
     * Plain Euclidean algorithm which fails if intermediate polynomials are not divisible at some step
     *
     * @param a poly
     * @param b poly
     * @return a list of polynomial remainders where the last element is GCD
     * @throws IllegalArgumentException if at some step intermediate polynomials are not divisible
     */
    public static PolynomialRemainders Euclid(final MutablePolynomial a,
                                              final MutablePolynomial b) {
        if (a.degree < b.degree)
            return Euclid(b, a);

        ArrayList<MutablePolynomial> prs = new ArrayList<>();
        prs.add(a.clone()); prs.add(b.clone());

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(prs);

        MutablePolynomial x = a, y = b, r;
        while (true) {
            //TODO replace with remainder!!!
            MutablePolynomial[] tmp = divideAndRemainder(x, y, true);
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
    public static PolynomialRemainders Euclid(final MutablePolynomial a,
                                              final MutablePolynomial b,
                                              final long modulus) {
        if (a.degree < b.degree)
            return Euclid(b, a, modulus);

        ArrayList<MutablePolynomial> prs = new ArrayList<>();
        prs.add(a.clone().modulus(modulus)); prs.add(b.clone().modulus(modulus));

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(prs);

        MutablePolynomial x = a, y = b, r;
        while (true) {
            //TODO replace with remainder!!!
            MutablePolynomial[] tmp = divideAndRemainder(x, y, modulus, true);
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
    public static PolynomialRemainders PolynomialEuclid(final MutablePolynomial a,
                                                        final MutablePolynomial b,
                                                        boolean primitivePRS) {
        if (a.degree < b.degree)
            return PolynomialEuclid(b, a, primitivePRS);


        if (a.isZero() || b.isZero()) return new PolynomialRemainders(a.clone(), b.clone());

        long aContent = a.content(), bContent = b.content();
        long contentGCD = gcd(aContent, bContent);
        MutablePolynomial aPP = a.clone().divide(aContent), bPP = b.clone().divide(bContent);

        ArrayList<MutablePolynomial> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        MutablePolynomial x = aPP, y = bPP, r;
        while (true) {
            //TODO replace with remainder!!!
            MutablePolynomial[] tmp = pseudoDivideAndRemainder(x, y, true);
            assert tmp != null;
            r = tmp[1];
            if (r.isZero())
                break;
            if (primitivePRS)
                r = r.primitivePart();
            prs.add(r);
            x = y;
            y = r;
        }
        PolynomialRemainders res = new PolynomialRemainders(prs);
        res.gcd().primitivePart().multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials which produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence where the last element is GCD
     */
    public static PolynomialRemainders SubresultantEuclid(final MutablePolynomial a,
                                                          final MutablePolynomial b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(a.clone(), b.clone());


        long aContent = a.content(), bContent = b.content();
        long contentGCD = gcd(aContent, bContent);
        MutablePolynomial aPP = a.clone().divide(aContent), bPP = b.clone().divide(bContent);

        ArrayList<MutablePolynomial> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            MutablePolynomial curr = prs.get(i);
            MutablePolynomial next = prs.get(i + 1);
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

            //TODO replace with remainder!!!
            MutablePolynomial q = pseudoDivideAndRemainder(curr, next, true)[1];
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
     * Representation for polynomial remainders sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders {
        public final ArrayList<MutablePolynomial> remainders;

        public PolynomialRemainders(MutablePolynomial... remainders) {
            this(new ArrayList<>(Arrays.asList(remainders)));
        }

        public PolynomialRemainders(ArrayList<MutablePolynomial> remainders) {
            this.remainders = remainders;
        }

        public MutablePolynomial gcd() {
            if (remainders.size() == 2 && remainders.get(1).isZero())
                return remainders.get(0);
            return remainders.get(remainders.size() - 1);
        }
    }

    /**
     * Modular GCD algorithm for polynomials
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static MutablePolynomial ModularGCD(MutablePolynomial a, MutablePolynomial b) {
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        long aContent = a.content(), bContent = b.content();
        long contentGCD = gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return MutablePolynomial.create(contentGCD);

        return ModularGCD0(a.clone().divide(aContent), b.clone().divide(bContent)).multiply(contentGCD);

    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static MutablePolynomial ModularGCD0(MutablePolynomial a, MutablePolynomial b) {
        if (a.degree < b.degree)
            return ModularGCD(b, a);

        long lcGCD = gcd(a.lc(), b.lc());
        double bound = Math.sqrt(a.degree + 1) * (1L << a.degree) * Math.max(a.norm(), b.norm()) * lcGCD;

        MutablePolynomial previousBase, base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            MutablePolynomial modularGCD = Euclid(a, b, prime).gcd();
            //clone if necessary
            if (modularGCD == a || modularGCD == b)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return MutablePolynomial.one();

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

            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            if ((double) basePrime >= 2 * bound || base.equals(previousBase)) {
                MutablePolynomial candidate = base.clone().symModulus(basePrime).primitivePart();
                //first check b since b is less degree
                MutablePolynomial[] div;
                div = divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
        }
    }

    /**
     * Computes GCD of two polynomials
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static MutablePolynomial PolynomialGCD(MutablePolynomial a, MutablePolynomial b) {
        return ModularGCD(a, b);
    }

    /**
     * Computes GCD of two polynomials modulo prime {@code modulus}
     *
     * @param a       the first polynomial
     * @param b       the second polynomial
     * @param modulus prime modulus
     * @return GCD of two polynomials
     */
    public static MutablePolynomial PolynomialGCD(MutablePolynomial a, MutablePolynomial b, long modulus) {
        return Euclid(a, b, modulus).gcd();
    }
}
