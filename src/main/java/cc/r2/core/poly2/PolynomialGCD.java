package cc.r2.core.poly2;


import cc.r2.core.number.primes.PrimesIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static cc.r2.core.poly2.DivisionWithRemainder.*;
import static cc.r2.core.poly2.LongArithmetics.*;

/**
 * Polynomial GCD and sub-resultant sequence for univariate polynomials with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class PolynomialGCD {
    private PolynomialGCD() {}

    /**
     * Euclidean algorithm
     *
     * @param a poly
     * @param b poly
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static <T extends IMutablePolynomial<T>> PolynomialRemainders<T> Euclid(final T a, final T b) {
        if (a.degree() < b.degree())
            return Euclid(b, a);

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
     * Runs extended Euclidean algorithm to compute {@code [gcd(a,b), x, y]} such that {@code x * a + y * b = gcd(a, b)}
     *
     * @param a the polynomial
     * @param b the polynomial
     * @return array of {@code [gcd(a,b), x, y]} such that {@code x * a + y * b = gcd(a, b)}
     */
    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T[] ExtendedEuclid(final T a, final T b) {
        T factory = a;
        T s = factory.createZero(), old_s = factory.createOne();
        T t = factory.createOne(), old_t = factory.createZero();
        T r = b, old_r = a;

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
        assert old_r.equals(a.clone().multiply(old_s).add(b.clone().multiply(old_t)));

        T[] result = a.arrayNewInstance(3);
        result[0] = old_r;
        result[1] = old_s;
        result[2] = old_t;
        return result;
    }

    /**
     * Euclidean algorithm for polynomials over Z that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to build primitive polynomial remainders or not
     * @return polynomial remainder sequence (the last element is GCD)
     */
    public static PolynomialRemainders<MutablePolynomialZ> PolynomialEuclid(final MutablePolynomialZ a,
                                                                            final MutablePolynomialZ b,
                                                                            boolean primitivePRS) {
        if (a.degree < b.degree)
            return PolynomialEuclid(b, a, primitivePRS);


        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());

        long aContent = a.content(), bContent = b.content();
        long contentGCD = LongArithmetics.gcd(aContent, bContent);
        MutablePolynomialZ aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<MutablePolynomialZ> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        MutablePolynomialZ x = aPP, y = bPP, r;
        while (true) {

            MutablePolynomialZ[] tmp = pseudoDivideAndRemainder(x, y, true);
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
        PolynomialRemainders<MutablePolynomialZ> res = new PolynomialRemainders<>(prs);
        res.gcd().primitivePart().multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials over Z that builds subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence (the last element is GCD)
     */
    public static PolynomialRemainders<MutablePolynomialZ> SubresultantEuclid(final MutablePolynomialZ a,
                                                                              final MutablePolynomialZ b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders<>(a.clone(), b.clone());


        long aContent = a.content(), bContent = b.content();
        long contentGCD = LongArithmetics.gcd(aContent, bContent);
        MutablePolynomialZ aPP = a.clone().divideOrNull(aContent), bPP = b.clone().divideOrNull(bContent);

        ArrayList<MutablePolynomialZ> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            MutablePolynomialZ curr = prs.get(i);
            MutablePolynomialZ next = prs.get(i + 1);
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

            MutablePolynomialZ q = pseudoDivideAndRemainder(curr, next, true)[1];
            if (q.isZero())
                break;

            q = q.divideOrNull(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders<MutablePolynomialZ> res = new PolynomialRemainders<>(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }


    /**
     * Polynomial remainder sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders<T extends IMutablePolynomial<T>> {
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
    }

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static MutablePolynomialZ ModularGCD(MutablePolynomialZ a, MutablePolynomialZ b) {
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        long aContent = a.content(), bContent = b.content();
        long contentGCD = gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return MutablePolynomialZ.create(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);

    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static MutablePolynomialZ ModularGCD0(MutablePolynomialZ a, MutablePolynomialZ b) {
        assert a.degree >= b.degree;

        long lcGCD = LongArithmetics.gcd(a.lc(), b.lc());
        double bound = Math.sqrt(a.degree + 1) * (1L << a.degree) * Math.max(a.norm2(), b.norm2()) * lcGCD;

        MutablePolynomialMod previousBase, base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            MutablePolynomialMod aMod = a.modulus(prime), bMod = b.modulus(prime);
            MutablePolynomialMod modularGCD = Euclid(aMod, bMod).gcd();
            //clone if necessary
            if (modularGCD == aMod || modularGCD == bMod)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return MutablePolynomialZ.one();

            //save the base
            if (base == null) {
                //make base monic and multiply lcGCD
                modularGCD.monic(lcGCD);
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
            long newBasePrime = safeMultiply(basePrime, prime);
            long monicFactor = modInverse(modularGCD.lc(), prime);
            long lcMod = mod(lcGCD, prime);
            for (int i = 0; i <= base.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                //long oth = mod(safeMultiply(mod(safeMultiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);

                long oth = modularGCD.multiplyMod(modularGCD.multiplyMod(modularGCD.data[i], monicFactor), lcMod);
                base.data[i] = ChineseRemainders(basePrime, prime, base.data[i], oth);
            }
            base = base.setModulusUnsafe(newBasePrime);
            basePrime = newBasePrime;

            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            if ((double) basePrime >= 2 * bound || base.equals(previousBase)) {
                MutablePolynomialZ candidate = base.normalSymmetricForm().primitivePart();
                //first check b since b is less degree
                MutablePolynomialZ[] div;
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
     * Returns GCD of two polynomials. Modular GCD algorithm is used for Z[x] and plain Euclid is used for Zp[x].
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("unchecked")
    public static <T extends MutablePolynomialAbstract<T>> T PolynomialGCD(T a, T b) {
        if (a instanceof MutablePolynomialZ)
            return (T) ModularGCD((MutablePolynomialZ) a, (MutablePolynomialZ) b);
        else
            return Euclid(a, b).gcd();
    }
}
