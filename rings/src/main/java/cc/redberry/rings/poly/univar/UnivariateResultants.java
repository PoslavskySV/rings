package cc.redberry.rings.poly.univar;

import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.Util;
import cc.redberry.rings.primes.PrimesIterator;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import static cc.redberry.rings.poly.univar.UnivariatePolynomial.asOverZp64;


/**
 * Various algorithms to compute (sub)resultants via Euclidean algorithm. Implementation is based on Gathen & Lücking,
 * "Subresultants revisited", https://doi.org/10.1016/S0304-3975(02)00639-4
 */
public final class UnivariateResultants {
    private UnivariateResultants() {}

    /**
     * Computes discriminant of polynomial and returns the result as a constant poly
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly DiscriminantAsPoly(Poly a) {
        if (a instanceof UnivariatePolynomialZp64)
            return (Poly) ((UnivariatePolynomialZp64) a).createConstant(Discriminant((UnivariatePolynomialZp64) a));
        else
            return (Poly) ((UnivariatePolynomial) a).createConstant(Discriminant((UnivariatePolynomial) a));
    }

    /**
     * Computes discriminant of polynomial
     */
    @SuppressWarnings("unchecked")
    public static <E> E Discriminant(UnivariatePolynomial<E> a) {
        Ring<E> ring = a.ring;
        E disc = ring.divideExact(Resultant(a, a.derivative()), a.lc());
        return ((a.degree * (a.degree - 1) / 2) % 2 == 1) ? ring.negate(disc) : disc;
    }

    /**
     * Computes discriminant of polynomial
     */
    @SuppressWarnings("unchecked")
    public static long Discriminant(UnivariatePolynomialZp64 a) {
        IntegersZp64 ring = a.ring;
        long disc = ring.divide(Resultant(a, a.derivative()), a.lc());
        return ((a.degree * (a.degree - 1) / 2) % 2 == 1) ? ring.negate(disc) : disc;
    }

    /**
     * Computes resultant of two polynomials and returns the result as a constant poly
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly ResultantAsPoly(Poly a, Poly b) {
        if (a instanceof UnivariatePolynomialZp64)
            return (Poly) ((UnivariatePolynomialZp64) a).createConstant(Resultant((UnivariatePolynomialZp64) a, (UnivariatePolynomialZp64) b));
        else
            return (Poly) ((UnivariatePolynomial) a).createConstant(Resultant((UnivariatePolynomial) a, (UnivariatePolynomial) b));
    }

    /**
     * Computes resultant of two polynomials
     */
    @SuppressWarnings("unchecked")
    public static <E> E Resultant(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        if (a.isOverField())
            return ClassicalPRS(a, b).resultant();
        else if (Util.isOverRationals(a))
            return (E) ResultantInQ((UnivariatePolynomial) a, (UnivariatePolynomial) b);

            //
            // Modular algorithm is really slow in practice (given no
            // any special bound) compared to Subresultant/Reduced
            //
            // else if (a.isOverZ())
            //     return (E) ModularResultant((UnivariatePolynomial) a, (UnivariatePolynomial) b);

        else
            return PrimitiveResultant(a, b, (p, q) -> SubresultantPRS(p, q).resultant());
    }

    /**
     * Computes resultant of two polynomials
     */
    public static long Resultant(UnivariatePolynomialZp64 a, UnivariatePolynomialZp64 b) {
        return ClassicalPRS(a, b).resultant();
    }

    private static <E> Rational<E> ResultantInQ(UnivariatePolynomial<Rational<E>> a,
                                                UnivariatePolynomial<Rational<E>> b) {
        Util.Tuple2<UnivariatePolynomial<E>, E>
                aZ = Util.toCommonDenominator(a),
                bZ = Util.toCommonDenominator(b);

        Ring<E> ring = aZ._1.ring;

        E resultant = Resultant(aZ._1, bZ._1);
        E den = ring.multiply(
                ring.pow(aZ._2, b.degree),
                ring.pow(bZ._2, a.degree));
        return new Rational<>(ring, resultant, den);
    }

    /**
     * Computes sequence of scalar subresultants.
     */
    public static <E> List<E> Subresultants(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        if (a.isOverField())
            return ClassicalPRS(a, b).getSubresultants();
        else
            return SubresultantPRS(a, b).getSubresultants();
    }

    static <E> E PrimitiveResultant(UnivariatePolynomial<E> a,
                                    UnivariatePolynomial<E> b,
                                    BiFunction<UnivariatePolynomial<E>, UnivariatePolynomial<E>, E> algorithm) {
        E ac = a.content(), bc = b.content();
        a = a.clone().divideExact(ac);
        b = b.clone().divideExact(bc);
        E r = algorithm.apply(a, b);
        Ring<E> ring = a.ring;
        r = ring.multiply(r, ring.pow(ac, b.degree));
        r = ring.multiply(r, ring.pow(bc, a.degree));
        return r;
    }

    /**
     * Modular algorithm for computing resultants over Z
     */
    public static BigInteger ModularResultant(UnivariatePolynomial<BigInteger> a,
                                              UnivariatePolynomial<BigInteger> b) {
        return PrimitiveResultant(a, b, UnivariateResultants::ModularResultant0);
    }

    private static BigInteger ModularResultant0(UnivariatePolynomial<BigInteger> a,
                                                UnivariatePolynomial<BigInteger> b) {
        // bound on the value of resultant
        BigInteger bound =
                UnivariatePolynomial.norm2(a).pow(b.degree)
                        .multiply(UnivariatePolynomial.norm2(b).pow(a.degree))
                        .shiftLeft(1);

        // aggregated CRT modulus
        BigInteger bModulus = null;
        BigInteger resultant = null;
        PrimesIterator primes = new PrimesIterator(1L << 15);
        while (true) {
            long prime = primes.take();
            BigInteger bPrime = BigInteger.valueOf(prime);
            IntegersZp zpRing = Rings.Zp(prime);
            UnivariatePolynomialZp64
                    aMod = asOverZp64(a.setRing(zpRing)),
                    bMod = asOverZp64(b.setRing(zpRing));

            if (aMod.degree != a.degree || bMod.degree != b.degree)
                continue;// unlucky prime

            long resultantMod = ClassicalPRS(aMod, bMod).resultant();
            BigInteger bResultantMod = BigInteger.valueOf(resultantMod);
            if (bModulus == null) {
                bModulus = bPrime;
                resultant = bResultantMod;
                continue;
            }

            if (!resultant.isZero() && resultantMod == 0)
                continue;// unlucky prime

            resultant = ChineseRemainders.ChineseRemainders(bModulus, bPrime, resultant, bResultantMod);
            bModulus = bModulus.multiply(BigInteger.valueOf(prime));

            if (bModulus.compareTo(bound) > 0)
                return Rings.Zp(bModulus).symmetricForm(resultant);
        }
    }

    /** Computes polynomial remainder sequence using classical division algorithm */
    public static PolynomialRemainderSequenceZp64 ClassicalPRS(UnivariatePolynomialZp64 a, UnivariatePolynomialZp64 b) {
        return new PolynomialRemainderSequenceZp64(a, b).run();
    }

    /** Computes polynomial remainder sequence using classical division algorithm */
    public static <E> PolynomialRemainderSequence<E> ClassicalPRS(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new ClassicalPolynomialRemainderSequence<>(a, b).run();
    }

    /** Computes polynomial remainder sequence using pseudo division algorithm */
    public static <E> PolynomialRemainderSequence<E> PseudoPRS(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new PseudoPolynomialRemainderSequence<>(a, b).run();
    }

    /** Computes polynomial remainder sequence using primitive division algorithm */
    public static <E> PolynomialRemainderSequence<E> PrimitivePRS(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new PrimitivePolynomialRemainderSequence<>(a, b).run();
    }

    /** Computes polynomial remainder sequence using reduced division algorithm */
    public static <E> PolynomialRemainderSequence<E> ReducedPRS(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new ReducedPolynomialRemainderSequence<>(a, b).run();
    }

    /** Computes subresultant polynomial remainder sequence */
    public static <E> PolynomialRemainderSequence<E> SubresultantPRS(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new SubresultantPolynomialRemainderSequence<>(a, b).run();
    }

    /**
     * Polynomial remainder sequence (PRS).
     */
    public static class APolynomialRemainderSequence<Poly extends IUnivariatePolynomial<Poly>> {
        /** Polynomial remainder sequence */
        public final List<Poly> remainders = new ArrayList<>();
        /** Quotients arised in PRS */
        public final List<Poly> quotients = new ArrayList<>();
        /** Initial polynomials */
        public final Poly a, b;

        public APolynomialRemainderSequence(Poly a, Poly b) {
            this.a = a;
            this.b = b;
        }

        /** The last element in PRS, that is the GCD */
        public final Poly lastRemainder() {
            return remainders.get(remainders.size() - 1);
        }

        /** The last element in PRS, that is the GCD */
        @SuppressWarnings("unchecked")
        public final Poly gcd() {
            if (a.isZero()) return b;
            if (b.isZero()) return a;

            if (a.isOverField())
                return lastRemainder().clone().monic();

            Poly r = lastRemainder().clone().primitivePartSameSign();
            return UnivariateGCD.PolynomialGCD(a.contentAsPoly(), b.contentAsPoly()).multiply(r);
        }

        public final int size() {
            return remainders.size();
        }
    }

    /**
     * Polynomial remainder sequence (PRS). It also implements abstract division rule, used to build PRS. At each step
     * of Euclidean algorithm the polynomials {@code qout, rem} and coefficients {@code alpha, beta} are computed so
     * that {@code alpha_i r_(i - 2) = quot_(i - 1) * r_(i - 1) + beta_i * r_i} where {@code {r_i} } is PRS.
     */
    public static abstract class PolynomialRemainderSequence<E> extends APolynomialRemainderSequence<UnivariatePolynomial<E>> {
        /** alpha coefficients */
        public final List<E> alphas = new ArrayList<>();
        /** beta coefficients */
        public final List<E> betas = new ArrayList<>();
        /** the ring */
        final Ring<E> ring;
        /** whether the first poly had smaller degree than the second */
        final boolean swap;

        PolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
            super(a, b);
            this.ring = a.ring;
            if (a.degree >= b.degree) {
                remainders.add(a);
                remainders.add(b);
                swap = false;
            } else {
                remainders.add(b);
                remainders.add(a);
                swap = a.degree % 2 == 1 && b.degree % 2 == 1; // both degrees are odd => odd permutation of Sylvester matrix
            }
        }

        /** compute alpha based on obtained so far PRS */
        abstract E nextAlpha();

        /** compute beta based on obtained so far PRS and newly computed remainder */
        abstract E nextBeta(UnivariatePolynomial<E> remainder);

        /** A single step of the Euclidean algorithm */
        private UnivariatePolynomial<E> step() {
            int i = remainders.size();
            UnivariatePolynomial<E>
                    dividend = remainders.get(i - 2).clone(),
                    divider = remainders.get(i - 1);

            E alpha = nextAlpha();
            dividend = dividend.multiply(alpha);

            UnivariatePolynomial<E>[] qd = UnivariateDivision.divideAndRemainder(dividend, divider, false);
            if (qd == null)
                throw new RuntimeException("exact division is not possible");

            UnivariatePolynomial<E> quotient = qd[0], remainder = qd[1];
            if (remainder.isZero())
                // remainder is zero => termination of the algorithm
                return remainder;

            E beta = nextBeta(remainder);
            remainder = remainder.divideExact(beta);

            alphas.add(alpha);
            betas.add(beta);
            quotients.add(quotient);
            remainders.add(remainder);
            return remainder;
        }

        /** Run all steps. */
        final PolynomialRemainderSequence<E> run() {
            if (lastRemainder().isZero())
                // on of the factors is zero
                return this;
            while (!step().isZero()) ; return this;
        }

        /** n_i - n_{i+1} */
        final int degreeDiff(int i) {
            return remainders.get(i).degree - remainders.get(i + 1).degree;
        }

        /** scalar subresultants */
        private final ArrayList<E> subresultants = new ArrayList<>();

        synchronized final void computeSubresultants() {
            if (!subresultants.isEmpty())
                return;

            List<E> subresultants = nonZeroSubresultants();
            if (swap) subresultants.replaceAll(ring::negate);
            this.subresultants.ensureCapacity(remainders.get(1).degree);
            for (int i = 0; i <= remainders.get(1).degree; ++i)
                this.subresultants.add(ring.getZero());
            for (int i = 1; i < remainders.size(); i++)
                this.subresultants.set(remainders.get(i).degree, subresultants.get(i - 1));
        }

        /** general setting for Fundamental Theorem of Resultant Theory */
        List<E> nonZeroSubresultants() {
            List<E> subresultants = new ArrayList<>();
            // largest subresultant
            E subresultant = ring.pow(remainders.get(1).lc(), degreeDiff(0));
            subresultants.add(subresultant);

            for (int i = 1; i < (remainders.size() - 1); ++i) {
                // computing (i+1)-th degree subresultant

                int di = degreeDiff(i);
                E rho = ring.pow(ring.multiply(remainders.get(i + 1).lc(), remainders.get(i).lc()), di);
                E den = ring.getOne();
                for (int j = 1; j <= i; ++j) {
                    rho = ring.multiply(rho, ring.pow(betas.get(j - 1), di));
                    den = ring.multiply(den, ring.pow(alphas.get(j - 1), di));
                }
                if ((di % 2) == 1 && (remainders.get(0).degree - remainders.get(i + 1).degree + i + 1) % 2 == 1)
                    rho = ring.negate(rho);
                subresultant = ring.multiply(subresultant, rho);
                subresultant = ring.divideExact(subresultant, den);

                subresultants.add(subresultant);
            }
            return subresultants;
        }

        /**
         * Gives a list of scalar subresultant where i-th list element is i-th subresultant.
         */
        public final List<E> getSubresultants() {
            if (subresultants.isEmpty())
                computeSubresultants();
            return subresultants;
        }

        /** Resultant of initial polynomials */
        public final E resultant() {
            return getSubresultants().get(0);
        }
    }

    /**
     * Classical division rule with alpha = beta = 1
     */
    private static final class ClassicalPolynomialRemainderSequence<E> extends PolynomialRemainderSequence<E> {
        ClassicalPolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
            super(a, b);
        }

        @Override
        final E nextAlpha() { return ring.getOne(); }

        @Override
        E nextBeta(UnivariatePolynomial<E> remainder) { return ring.getOne(); }

        @Override
        List<E> nonZeroSubresultants() {
            List<E> subresultants = new ArrayList<>();
            // largest subresultant
            E subresultant = ring.pow(remainders.get(1).lc(), degreeDiff(0));
            subresultants.add(subresultant);

            for (int i = 1; i < (remainders.size() - 1); ++i) {
                // computing (i+1)-th degree subresultant

                int di = degreeDiff(i);
                E rho = ring.pow(ring.multiply(remainders.get(i + 1).lc(), remainders.get(i).lc()), di);
                if ((di % 2) == 1 && (remainders.get(0).degree - remainders.get(i + 1).degree + i + 1) % 2 == 1)
                    rho = ring.negate(rho);
                subresultant = ring.multiply(subresultant, rho);
                subresultants.add(subresultant);
            }
            return subresultants;
        }
    }

    private static class PseudoPolynomialRemainderSequence<E> extends PolynomialRemainderSequence<E> {
        PseudoPolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
            super(a, b);
        }

        @Override
        final E nextAlpha() {
            int i = remainders.size();
            E lc = remainders.get(i - 1).lc();
            int deg = remainders.get(i - 2).degree - remainders.get(i - 1).degree;
            return ring.pow(lc, deg + 1);
        }

        @Override
        E nextBeta(UnivariatePolynomial<E> remainder) {
            return ring.getOne();
        }
    }

    /**
     * Reduced pseudo-division
     */
    private static final class ReducedPolynomialRemainderSequence<E> extends PseudoPolynomialRemainderSequence<E> {
        ReducedPolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
            super(a, b);
        }

        @Override
        E nextBeta(UnivariatePolynomial<E> remainder) {
            return alphas.isEmpty() ? ring.getOne() : alphas.get(alphas.size() - 1);
        }

        @Override
        List<E> nonZeroSubresultants() {
            List<E> subresultants = new ArrayList<>();
            // largest subresultant
            E subresultant = ring.pow(remainders.get(1).lc(), degreeDiff(0));
            subresultants.add(subresultant);

            for (int i = 1; i < (remainders.size() - 1); ++i) {
                // computing (i+1)-th degree subresultant
                int di = degreeDiff(i);
                E rho = ring.pow(remainders.get(i + 1).lc(), di);
                E den = ring.pow(remainders.get(i).lc(), degreeDiff(i - 1) * di);
                subresultant = ring.multiply(subresultant, rho);
                subresultant = ring.divideExact(subresultant, den);
                if ((di % 2) == 1 && (remainders.get(0).degree - remainders.get(i + 1).degree + i + 1) % 2 == 1)
                    subresultant = ring.negate(subresultant);
                subresultants.add(subresultant);
            }
            return subresultants;
        }
    }

    /**
     * Primitive pseudo-division
     */
    private static final class PrimitivePolynomialRemainderSequence<E> extends PseudoPolynomialRemainderSequence<E> {
        PrimitivePolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) { super(a, b); }

        @Override
        E nextBeta(UnivariatePolynomial<E> remainder) { return remainder.content(); }
    }

    /**
     * Subresultant sequence
     */
    private static final class SubresultantPolynomialRemainderSequence<E> extends PseudoPolynomialRemainderSequence<E> {
        final List<E> psis = new ArrayList<>();

        SubresultantPolynomialRemainderSequence(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) { super(a, b); }

        @Override
        E nextBeta(UnivariatePolynomial<E> remainder) {
            int i = remainders.size();
            UnivariatePolynomial<E> prem = remainders.get(i - 2);
            E lc = i == 2 ? ring.getOne() : prem.lc();
            E psi;
            if (i == 2)
                psi = ring.getNegativeOne();
            else {
                E prevPsi = psis.get(psis.size() - 1);
                int deg = remainders.get(i - 3).degree - remainders.get(i - 2).degree;
                E f = ring.pow(ring.negate(lc), deg);
                if (1 - deg < 0)
                    psi = ring.divideExact(f, ring.pow(prevPsi, deg - 1));
                else
                    psi = ring.multiply(f, ring.pow(prevPsi, 1 - deg));
            }
            psis.add(psi);
            return ring.negate(ring.multiply(lc, ring.pow(psi, remainders.get(i - 2).degree - remainders.get(i - 1).degree)));
        }

        private int eij(int i, int j) {
            int e = degreeDiff(j - 1);
            for (int k = j; k <= i; ++k)
                e *= 1 - degreeDiff(k);
            return e;
        }

        @Override
        List<E> nonZeroSubresultants() {
            List<E> subresultants = new ArrayList<>();
            // largest subresultant
            E subresultant = ring.pow(remainders.get(1).lc(), degreeDiff(0));
            subresultants.add(subresultant);

            for (int i = 1; i < (remainders.size() - 1); ++i) {
                // computing (i+1)-th degree subresultant

                int di = degreeDiff(i);
                E rho = ring.pow(remainders.get(i + 1).lc(), di);
                E den = ring.getOne();
                for (int k = 1; k <= i; ++k) {
                    int deg = -di * eij(i - 1, k);
                    if (deg >= 0)
                        rho = ring.multiply(rho, ring.pow(remainders.get(k).lc(), deg));
                    else
                        den = ring.multiply(den, ring.pow(remainders.get(k).lc(), -deg));
                }
                subresultant = ring.multiply(subresultant, rho);
                subresultant = ring.divideExact(subresultant, den);
                subresultants.add(subresultant);
            }
            return subresultants;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Classical division rule for polynomials over Zp
     */
    public static final class PolynomialRemainderSequenceZp64 extends APolynomialRemainderSequence<UnivariatePolynomialZp64> {
        /** the ring */
        final IntegersZp64 ring;
        /** whether the first poly had smaller degree than the second */
        final boolean swap;

        PolynomialRemainderSequenceZp64(UnivariatePolynomialZp64 a, UnivariatePolynomialZp64 b) {
            super(a, b);
            this.ring = a.ring;
            if (a.degree >= b.degree) {
                remainders.add(a);
                remainders.add(b);
                swap = false;
            } else {
                remainders.add(b);
                remainders.add(a);
                swap = a.degree % 2 == 1 && b.degree % 2 == 1; // both degrees are odd => odd permutation of Sylvester matrix
            }
        }

        /** A single step of the Euclidean algorithm */
        private UnivariatePolynomialZp64 step() {
            int i = remainders.size();
            UnivariatePolynomialZp64
                    dividend = remainders.get(i - 2).clone(),
                    divider = remainders.get(i - 1);

            UnivariatePolynomialZp64[] qd = UnivariateDivision.divideAndRemainder(dividend, divider, false);
            if (qd == null)
                throw new RuntimeException("exact division is not possible");

            UnivariatePolynomialZp64 quotient = qd[0], remainder = qd[1];
            if (remainder.isZero())
                // remainder is zero => termination of the algorithm
                return remainder;

            quotients.add(quotient);
            remainders.add(remainder);
            return remainder;
        }

        /** Run all steps. */
        private PolynomialRemainderSequenceZp64 run() { while (!step().isZero()) ; return this;}

        /** n_i - n_{i+1} */
        final int degreeDiff(int i) {
            return remainders.get(i).degree - remainders.get(i + 1).degree;
        }

        /** scalar subresultants */
        private final TLongArrayList subresultants = new TLongArrayList();

        synchronized final void computeSubresultants() {
            if (!subresultants.isEmpty())
                return;

            TLongArrayList subresultants = nonZeroSubresultants();
            if (swap)
                for (int i = 0; i < subresultants.size(); ++i)
                    subresultants.set(i, ring.negate(subresultants.get(i)));

            this.subresultants.ensureCapacity(remainders.get(1).degree);
            for (int i = 0; i <= remainders.get(1).degree; ++i)
                this.subresultants.add(0L);
            for (int i = 1; i < remainders.size(); i++)
                this.subresultants.set(remainders.get(i).degree, subresultants.get(i - 1));
        }

        /** general setting for Fundamental Theorem of Resultant Theory */
        TLongArrayList nonZeroSubresultants() {
            TLongArrayList subresultants = new TLongArrayList();
            // largest subresultant
            long subresultant = ring.powMod(remainders.get(1).lc(), degreeDiff(0));
            subresultants.add(subresultant);

            for (int i = 1; i < (remainders.size() - 1); ++i) {
                // computing (i+1)-th degree subresultant

                int di = degreeDiff(i);
                long rho = ring.powMod(ring.multiply(remainders.get(i + 1).lc(), remainders.get(i).lc()), di);
                if ((di % 2) == 1 && (remainders.get(0).degree - remainders.get(i + 1).degree + i + 1) % 2 == 1)
                    rho = ring.negate(rho);
                subresultant = ring.multiply(subresultant, rho);
                subresultants.add(subresultant);
            }
            return subresultants;
        }

        /**
         * Gives a list of scalar subresultant where i-th list element is i-th subresultant.
         */
        public final TLongArrayList getSubresultants() {
            if (subresultants.isEmpty())
                computeSubresultants();
            return subresultants;
        }

        /** Resultant of initial polynomials */
        public final long resultant() {
            return getSubresultants().get(0);
        }
    }
}