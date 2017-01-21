package cc.r2.core.polynomial;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.LongArithmetics.*;
import static cc.r2.core.polynomial.MutablePolynomial.multiplyModClassical;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 14/01/2017.
 */
@Ignore
public class FFT {

    static int revbin(int in, int bits) {
        return Integer.reverse(in) >>> (32 - bits);
    }

    static void revbin(long[] target, final int bits, final long[] data) {
        for (int i = 0; i < data.length; i++)
            target[revbin(i, bits)] = data[i];
    }

    static long[] FastFourier(final int bits, final long[] data, final long[] rootPowers, final long modulus) {
        long[] rev = new long[data.length];
        revbin(rev, bits, data);
        return FastFourier0(bits, rev, rootPowers, modulus);
    }

    static long[] FastFourier0(int bits, long[] revbinData, long[] rootPowers, long modulus) {
        int n = revbinData.length;
        for (int s = 1; s <= bits; ++s) {
            int m = 1 << s;
            for (int k = 0; k <= n - 1; k += m) {
                for (int j = 0; j <= m / 2 - 1; ++j) {
                    long t = multiplyMod(rootPowers[j * n / m], revbinData[k + j + m / 2], modulus);
                    long u = revbinData[k + j];
                    revbinData[k + j] = addMod(u, t, modulus);
                    revbinData[k + j + m / 2] = addMod(u, -t, modulus);
                }
            }
        }
        return revbinData;
    }

    public static final class PolynomialFFT {
        final long[] fft;
        final long[] fftInverse;

        final long modulus;
        final long[] rootPowers;
        final long[] rootPowersInv;
        final int bits;
        int degree;
        final long nInversed;

        public PolynomialFFT(int bits, long[] rootPowers, long[] rootPowersInv, int degree, long[] coefficients, long modulus) {
            assert rootPowers.length == 1 << bits;

            this.bits = bits;
            this.degree = degree;
            this.modulus = modulus;
            this.rootPowers = rootPowers;
            this.rootPowersInv = rootPowersInv;
            this.nInversed = modInverse(rootPowers.length, modulus);

            fft = new long[rootPowers.length];
            revbin(fft, bits, Arrays.copyOf(coefficients, rootPowers.length));
            fftInverse = fft.clone();
            FastFourier0(bits, fft, rootPowers, modulus);
            FastFourier0(bits, fftInverse, rootPowersInv, modulus);
        }

        private PolynomialFFT(long[] fft, long[] fftInverse, long modulus, long[] rootPowers, long[] rootPowersInv, int bits, int degree, long nInversed) {
            this.fft = fft;
            this.fftInverse = fftInverse;
            this.modulus = modulus;
            this.rootPowers = rootPowers;
            this.rootPowersInv = rootPowersInv;
            this.bits = bits;
            this.degree = degree;
            this.nInversed = nInversed;
        }

        private void checkCompatibility(PolynomialFFT other) {
            if (bits != other.bits || modulus != other.modulus || rootPowers[0] != other.rootPowers[0])
                throw new IllegalArgumentException();
        }

        PolynomialFFT createOne() {
            final long[] fft = new long[this.fft.length];
            final long[] fftInverse = new long[this.fft.length];
            Arrays.fill(fft, 1);
            Arrays.fill(fftInverse, 1);
            return new PolynomialFFT(fft, fftInverse, modulus, rootPowers, rootPowersInv, bits, 0, nInversed);
        }

        PolynomialFFT createZero() {
            final long[] fft = new long[this.fft.length];
            final long[] fftInverse = new long[this.fft.length];
            Arrays.fill(fft, 0);
            Arrays.fill(fftInverse, 0);
            return new PolynomialFFT(fft, fftInverse, modulus, rootPowers, rootPowersInv, bits, 0, nInversed);
        }

        PolynomialFFT createMonomial(int degree) {
            final long[] fft = new long[this.fft.length];
            final long[] fftInverse = new long[this.fft.length];
            for (int i = 0; i < fft.length; i++) {
                int rootPower = (degree * i) % rootPowers.length;
                fft[i] = rootPowers[rootPower];
                fftInverse[i] = rootPowersInv[rootPower];
            }
            return new PolynomialFFT(fft, fftInverse, modulus, rootPowers, rootPowersInv, bits, 0, nInversed);
        }

        PolynomialFFT add(PolynomialFFT other) {
            checkCompatibility(other);
            for (int i = 0; i < fft.length; i++) {
                fft[i] = mod(fft[i] + other.fft[i], modulus);
                fftInverse[i] = mod(fftInverse[i] + other.fftInverse[i], modulus);
            }
            degree = Math.max(degree, other.degree);
            return this;
        }

        PolynomialFFT subtract(PolynomialFFT other) {
            checkCompatibility(other);
            for (int i = 0; i < fft.length; ++i) {
                fft[i] = mod(fft[i] - other.fft[i], modulus);
                fftInverse[i] = mod(fftInverse[i] - other.fftInverse[i], modulus);
            }
            degree = Math.max(degree, other.degree);
            return this;
        }

        PolynomialFFT multiply(PolynomialFFT other) {
            checkCompatibility(other);
            if (degree + other.degree > fft.length)
                throw new IllegalArgumentException("Not enough space: space = " + (fft.length) + " this.degree = " + degree + " oth.degree = " + other.degree);
            for (int i = 0; i < fft.length; i++) {
                fft[i] = mod(fft[i] * other.fft[i], modulus);
                fftInverse[i] = mod(fftInverse[i] * other.fftInverse[i], modulus);
            }
            degree += other.degree;
            return this;
        }

        PolynomialFFT multiply(long val) {
            for (int i = 0; i < fft.length; i++) {
                fft[i] = mod(fft[i] * val, modulus);
                fftInverse[i] = mod(fftInverse[i] * val, modulus);
            }
            return this;
        }

        PolynomialFFT square() {
            checkCompatibility(this);
            if (degree + degree > fft.length)
                throw new IllegalArgumentException("Not enough space");
            for (int i = 0; i < fft.length; i++) {
                fft[i] = mod(fft[i] * fft[i], modulus);
                fftInverse[i] = mod(fftInverse[i] * fftInverse[i], modulus);
            }
            degree += degree;
            return this;
        }

        PolynomialFFT reverse() {
            for (int i = 0; i < fft.length; ++i) {
                long tmp = fft[i];
                int rootPower = (i * degree) % rootPowers.length;
                fft[i] = mod(fftInverse[i] * rootPowers[rootPower], modulus);
                fftInverse[i] = mod(tmp * rootPowersInv[rootPower], modulus);
            }
            return this;
        }

        PolynomialFFT cut(int newDegree) {
            if (newDegree >= degree)
                return this;

            long[] newCoefficients = toPoly().cut(newDegree).data;
//            long[] newCoefficients = new long[fft.length];
//            Arrays.fill(newCoefficients, 3);
            revbin(fft, bits, newCoefficients);
            System.arraycopy(fft, 0, fftInverse, 0, fft.length);
            FastFourier0(bits, fft, rootPowers, modulus);
            FastFourier0(bits, fftInverse, rootPowersInv, modulus);

            this.degree = newDegree;
            return this;
        }

        PolynomialFFT shiftRight(int n) {
            for (int i = 0; i < fft.length; i++) {
                int rootPower = (i * n) % rootPowers.length;
                fft[i] = mod(fft[i] * rootPowers[rootPower], modulus);
                fftInverse[i] = mod(fftInverse[i] * rootPowersInv[rootPower], modulus);
            }
            degree += n;
            return this;
        }

        private long[] polyCoefficients() {
            long[] coefficients = FastFourier(bits, fft, rootPowersInv, modulus);
            for (int i = 0; i < coefficients.length; i++)
                coefficients[i] = mod(coefficients[i] * nInversed, modulus);
            return coefficients;
        }

        MutablePolynomial toPoly() {
            return MutablePolynomial.create(polyCoefficients());
        }

        @Override
        public PolynomialFFT clone() {
            return new PolynomialFFT(fft.clone(), fftInverse.clone(), modulus, rootPowers, rootPowersInv, bits, degree, nInversed);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PolynomialFFT that = (PolynomialFFT) o;

            if (modulus != that.modulus) return false;
            if (bits != that.bits) return false;
            if (degree != that.degree) return false;
            if (nInversed != that.nInversed) return false;
            if (rootPowers[1] != that.rootPowers[1]) return false;
            if (!Arrays.equals(fft, that.fft)) return false;
            if (!Arrays.equals(fftInverse, that.fftInverse)) return false;
            return Arrays.equals(rootPowersInv, that.rootPowersInv);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(fft);
            result = 31 * result + Arrays.hashCode(fftInverse);
            result = 31 * result + (int) (modulus^(modulus >>> 32));
            result = 31 * result + (int) (modulus^(rootPowers[1] >>> 32));
            result = 31 * result + bits;
            result = 31 * result + degree;
            result = 31 * result + (int) (nInversed^(nInversed >>> 32));
            return result;
        }
    }

    static PolynomialFFT toFFT(MutablePolynomial a) {
        return new PolynomialFFT(nPower, rootPowers, rootInvPowers, a.degree, a.data, pModulus);
    }

    static PolynomialFFT cut0(PolynomialFFT poly, int newDegree) {
        if (newDegree >= poly.degree)
            return poly;

        long[] newCoefficients = poly.toPoly().cut(newDegree).data;
        revbin(poly.fft, poly.bits, newCoefficients);
        System.arraycopy(poly.fft, 0, poly.fftInverse, 0, poly.fft.length);
        FastFourier0(poly.bits, poly.fft, poly.rootPowers, poly.modulus);
        FastFourier0(poly.bits, poly.fftInverse, poly.rootPowersInv, poly.modulus);

        poly.degree = newDegree;
        return poly;
    }

    static PolynomialFFT cut1(PolynomialFFT poly, int newDegree) {
        if (newDegree >= poly.degree)
            return poly;

        for (int i = 0; i < poly.fft.length; i++) {
            int rootPower = ((newDegree) * i) % poly.rootPowers.length;
            poly.fft[i] = mod(mod(poly.fft[i], poly.rootPowers[rootPower]), pModulus);
            poly.fftInverse[i] = mod(mod(poly.fftInverse[i], poly.rootPowersInv[rootPower]), pModulus);
        }
        poly.degree = newDegree;
        return poly;
    }

    @Test
    public void test_poly_cut() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 122, 3, 4, 5, 6, 6, 1, 2, 34);
        PolynomialFFT aFFt = toFFT(a);
        System.out.println(aFFt.clone().cut(5).toPoly());
        System.out.println(cut1(aFFt.clone(), 5).toPoly());
    }

    @Test
    public void test_poly1() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3);
        MutablePolynomial b = MutablePolynomial.create(1, 2, 3, 4, 5);

        PolynomialFFT fftA = toFFT(a);
        PolynomialFFT fftB = toFFT(b);

        Assert.assertEquals(a.clone().add(b, pModulus), fftA.clone().add(fftB).toPoly());

        Assert.assertEquals(a.clone().subtract(b, pModulus), fftA.clone().subtract(fftB).toPoly());
        Assert.assertEquals(toFFT(a.clone().subtract(b, pModulus)), fftA.clone().subtract(fftB));

        Assert.assertEquals(a.clone().multiply(b, pModulus), fftA.clone().multiply(fftB).toPoly());
        Assert.assertEquals(toFFT(a.clone().multiply(b, pModulus)), fftA.clone().multiply(fftB));
    }


    @Test
    public void test_poly2() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5), rev = a.clone().reverse();
        PolynomialFFT fft = toFFT(a);
        PolynomialFFT fftRev = toFFT(rev);
        Assert.assertEquals(fft.reverse(), fftRev);
    }

    @Test
    public void test_poly3() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5);
        PolynomialFFT fft = toFFT(a);
        for (int i = 0; i < 3; i++) {
            Assert.assertEquals(a.clone().cut(i), fft.clone().cut(i).toPoly());
            Assert.assertEquals(toFFT(a.clone().cut(i)), fft.clone().cut(i));
        }
    }

    @Test
    public void test_poly4() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5);
        PolynomialFFT fft = new PolynomialFFT(nPower, rootPowers, rootInvPowers, a.degree, a.data, pModulus);
        for (int i = 0; i < 10; i++) {
            Assert.assertEquals(a.clone().shiftRight(i), fft.clone().shiftRight(i).toPoly());
            Assert.assertEquals(toFFT(a.clone().shiftRight(i)), fft.clone().shiftRight(i));
        }
    }

    public static PolynomialFFT remainderMonomial(PolynomialFFT dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).cut(xDegree - 1);
    }

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomialFFT {
        final PolynomialFFT poly;

        private InverseModMonomialFFT(PolynomialFFT poly) {
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<PolynomialFFT> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public PolynomialFFT getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() > r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            PolynomialFFT gPrev = currentSize == 0 ? poly.createOne() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                PolynomialFFT tmp = gPrev.clone().multiply(2L).subtract(gPrev.clone().square().multiply(poly));
                inverses.add(gPrev = remainderMonomial(tmp, 1 << i, false));
            }
            return gPrev;
        }
    }

    public static InverseModMonomialFFT fastDivisionPreConditioning(PolynomialFFT divider) {
        return new InverseModMonomialFFT(divider.clone().reverse());
    }

    @Test
    public void test_poly6() throws Exception {
        MutablePolynomial b = MutablePolynomial.create(1, 2, 3, 4, 5, 1);
        PolynomialFFT fftB = toFFT(b);
        InverseModMonomialFFT inv = new InverseModMonomialFFT(fftB.clone().reverse());
        InverseModMonomial inv0 = new InverseModMonomial(b.clone().reverse(), pModulus);
        for (int i = 1; i < 33; i++) {
            assertEquals(inv0.getInverse(i), inv.getInverse(i).toPoly());
            assertInverseModMonomial(fftB.clone().reverse(), inv.getInverse(i), i, pModulus);
        }
    }

    static void assertInverseModMonomial(PolynomialFFT poly, PolynomialFFT invMod, int monomialDegree, long modulus) {
        Assert.assertTrue(PolynomialArithmetics.polyMod(poly.clone().multiply(invMod).toPoly(), MutablePolynomial.createMonomial(1, monomialDegree), modulus, true).isOne());
    }

    static PolynomialFFT[] divideAndRemainderFast0(PolynomialFFT dividend,
                                                   PolynomialFFT divider,
                                                   InverseModMonomialFFT invRevMod,
                                                   boolean copy) {
        if (dividend.degree < divider.degree)
            return new PolynomialFFT[]{dividend.createZero(), copy ? dividend.clone() : dividend};

        int m = dividend.degree - divider.degree;
        PolynomialFFT quot = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1)), m + 1, false).reverse();
        if (quot.degree < m)
            quot.shiftRight(m - quot.degree);
        quot.degree = dividend.degree - divider.degree;
        PolynomialFFT rem = (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(quot));
        rem.degree = Math.min(rem.degree, divider.degree - 1);
        return new PolynomialFFT[]{quot, rem};
    }

    static PolynomialFFT remainder(PolynomialFFT dividend,
                                   PolynomialFFT divider,
                                   InverseModMonomialFFT invRevMod,
                                   boolean copy) {
        return divideAndRemainderFast0(dividend, divider, invRevMod, copy)[1];
    }

    @Test
    public void test_poly5() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 67, 7);
        MutablePolynomial b = MutablePolynomial.create(1, 2, 3, 4, 1);

        PolynomialFFT fftA = new PolynomialFFT(nPower, rootPowers, rootInvPowers, a.degree, a.data, pModulus);
        PolynomialFFT fftB = new PolynomialFFT(nPower, rootPowers, rootInvPowers, b.degree, b.data, pModulus);

        InverseModMonomialFFT bInv = fastDivisionPreConditioning(fftB.clone());
        PolynomialFFT[] qd = divideAndRemainderFast0(fftA.clone(), fftB.clone(), bInv, true);
        MutablePolynomial[] qdPoly = {qd[0].toPoly(), qd[1].toPoly()};
        System.out.println(Arrays.toString(qdPoly));
        System.out.println(Arrays.toString(DivideAndRemainder.divideAndRemainder(a, b, pModulus, true)));
    }


    public static PolynomialFFT polyPowMod(final PolynomialFFT base, long exponent, PolynomialFFT polyModulus, InverseModMonomialFFT invMod, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return base.createOne();

        PolynomialFFT result = base.createOne();
        PolynomialFFT k2p = remainder(base, polyModulus, invMod, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = remainder(result.multiply(k2p), polyModulus, invMod, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = remainder(k2p.multiply(k2p), polyModulus, invMod, false);
        }
    }


    public static ArrayList<MutablePolynomial> ExponentsNaive(MutablePolynomial polyModulus, long modulus, int n) {
        MutablePolynomial exponent = MutablePolynomial.create(0, 1);
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        InverseModMonomial invMod = DivideAndRemainder.fastDivisionPreConditioning(polyModulus, modulus);
        for (int j = 0; j < n; j++) {
            exps.add(exponent);
            exponent = PolynomialArithmetics.polyPowMod(exponent, modulus, polyModulus, invMod, modulus, true);
        }
        return exps;
    }

    public static ArrayList<PolynomialFFT> ExponentsNaiveFFT(MutablePolynomial polyModulus, long modulus, int n) {
        PolynomialFFT poly = toFFT(polyModulus);
        InverseModMonomialFFT invMod = fastDivisionPreConditioning(poly);
        return ExponentsNaiveFFT(poly, invMod, modulus, n);
    }

    public static ArrayList<PolynomialFFT> ExponentsNaiveFFT(PolynomialFFT poly, InverseModMonomialFFT invMod, long modulus, int n) {
        PolynomialFFT exponent = toFFT(MutablePolynomial.create(0, 1));
        ArrayList<PolynomialFFT> exps = new ArrayList<>();
        for (int j = 0; j < n; j++) {
            exps.add(exponent);
            exponent = polyPowMod(exponent, modulus, poly, invMod, true);
        }
        return exps;
    }

//     -----
//             1096502685
//             5163977853
//             -----
//             1069429689
//             5093414936
//             -----
//             1064802460
//             5040595146


    @Test
    public void sad() throws Exception {
        MutablePolynomial a = RandomPolynomials.randomMonicPoly(500, pModulus, new Well1024a());
        MutablePolynomial b = RandomPolynomials.randomMonicPoly(500, pModulus, new Well1024a());
        MutablePolynomial c = RandomPolynomials.randomMonicPoly(100, pModulus, new Well1024a());

        PolynomialFFT dov = toFFT(c);
        MutablePolynomial e = remainder(toFFT(a).multiply(toFFT(b)), dov, fastDivisionPreConditioning(dov), false).toPoly();
        System.out.println(e);
        System.out.println(DivideAndRemainder.remainder(a.multiply(b, pModulus), c, pModulus, true));

    }

    @Test
    public void name() throws Exception {
        MutablePolynomial poly = RandomPolynomials.randomMonicPoly(250, pModulus, new Well1024a());

        int nExps = 100;
        for (int i = 0; i < 1000; i++) {
            System.out.println(" ----- ");
            long start = System.nanoTime();
            ArrayList<MutablePolynomial> a = ExponentsNaive(poly, pModulus, nExps);
            System.out.println(System.nanoTime() - start);

            PolynomialFFT poly0 = toFFT(poly);
            InverseModMonomialFFT invMod = fastDivisionPreConditioning(poly0);

            start = System.nanoTime();
            ArrayList<PolynomialFFT> b = ExponentsNaiveFFT(poly0, invMod, pModulus, nExps);
            System.out.println(System.nanoTime() - start);

//            for (int j = 0; j < nExps; j++) {
//                assertEquals(a.get(j), b.get(j));
//            }
        }


    }

    static long[] fft_3(int r, long[] as, long val, long modulus) {
        long[] w = new long[as.length];
        long[] rev = new long[as.length];
        for (int i = 0; i < as.length; i++) {
            w[i] = powMod(val, i, modulus);
            rev[revbin(i, r)] = as[i];
        }
        return fft_30(r, rev, w, modulus);
    }

    static long[] fft_3(final int r, final long[] as, final long[] rootPowers, final long modulus) {
        long[] rev = new long[as.length];
        for (int i = 0; i < as.length; i++) {
            rev[revbin(i, r)] = as[i];
        }
        return fft_30(r, rev, rootPowers, modulus);
    }

    static long[] fft_3_mod(int n, int r, long[] as, long[] w, long modulus) {
        for (int s = 1; s <= r; ++s) {
            int m = 1 << s;
            for (int k = 0; k <= n - 1; k += m) {
                long o = w[0];
                for (int j = 0; j <= m / 2 - 1; ++j) {
                    long t = multiplyMod(w[j * n / m], as[k + j + m / 2], modulus);
                    long u = as[k + j];
                    as[k + j] = addMod(u, t, modulus);
                    as[k + j + m / 2] = addMod(u, -t, modulus);
                }
            }
        }
        return as;
    }

    static int counter = 0;

    @Test
    public void asdasd() throws Exception {

        System.out.println(124214L % pModulus);
        System.out.println((124214L - pModulus) % pModulus);
    }


    static long[] fft_30(final int r, final long[] as, final long[] w, final long modulus) {
        int n = as.length;
        for (int s = 1; s <= r; ++s) {
            int m = 1 << s;
            int nm = n / m;
            for (int k = 0; k < n; k += m) {
                for (int j = k, jTo = k + m / 2; j < jTo; ++j) {
                    ++counter;
                    long t = (w[(j - k) * nm] * as[j + m / 2]);
//                    if (t > modulus)
//                        t -= modulus;
                    t = t % modulus;
                    long u = as[j];
                    as[j] = u + t;
                    as[j + m / 2] = u - t;
                }
            }
        }
        return as;
    }


    @Test
    public void asdasdasdasd() throws Exception {
        long[] a = RandomPolynomials.randomLongArray(500, 100, new Well1024a());
        long[] b = RandomPolynomials.randomLongArray(500, 100, new Well1024a());
        a = Arrays.copyOf(a, n);
        b = Arrays.copyOf(b, n);

        long[] a_fft = fft_3(nPower, a, rootPowers, pModulus);
        counter += n;
        long[] b_fft = fft_3(nPower, b, rootPowers, pModulus);
        counter += n;
        long[] res = new long[a.length];
        for (int i = 0; i < res.length; i++)
            res[i] = (a_fft[i] * b_fft[i]) % pModulus;
        counter += n;

        res = fft_3(nPower, res, rootInvPowers, pModulus);
        counter += n;
        for (int i = 0; i < res.length; i++)
            res[i] = multiplyMod(res[i], nInv, pModulus);
        counter += n;

        System.out.println(500 * 500);
        System.out.println(counter);
    }


    //
//    @Test
//    public void name() throws Exception {
//        System.out.println();
//    }

    //
//    static final long pModulus = 17;
//    static final long root = 3;
//    static final long rootInv = LongArithmetics.modInverse(root, pModulus);
//    static final int n = 16;
//    static final int nPower = 4;
//    static final long nInv = modInverse(n, pModulus);

//    static final long pModulus = 257;
//    static final long root = 9;
//    static final long rootInv = LongArithmetics.modInverse(root, pModulus);
//    static final int n = 128;
//    static final int nPower = 7;
//    static final long nInv = modInverse(n, pModulus);

    static final long pModulus = 1022977;
    static final long root = 78;
    static final long rootInv = LongArithmetics.modInverse(root, pModulus);
    static final int n = 1024;
    static final int nPower = 10;
    static final long nInv = modInverse(n, pModulus);


//    static final long pModulus = 2021377;
//    static final long root = 1233;
//    static final long rootInv = LongArithmetics.modInverse(root, pModulus);
//    static final int n = 1 << 11;
//    static final int nPower = 11;
//    static final long nInv = modInverse(n, pModulus);

//    static final long pModulus = 7340033;
//    static final long root = 5;
//    static final long rootInv = LongArithmetics.modInverse(root, pModulus);
//    static final int n = 1 << 20;
//    static final int nPower = 20;
//    static final long nInv = modInverse(n, pModulus);


    static final long[] rootPowers;
    static final long[] rootInvPowers;

    static {
        rootPowers = new long[n];
        rootInvPowers = new long[n];
        for (int i = 0; i < n; i++) {
            rootPowers[i] = powMod(root, i, pModulus);
            rootInvPowers[i] = powMod(rootInv, i, pModulus);
        }
    }

    static boolean isPRoot(long root, int exp, long modulus) {
        if (powMod(root, exp, modulus) != 1)
            return false;
        for (int i = 1; i < exp; i++)
            if (powMod(root, i, modulus) == 1)
                return false;
        return true;
    }

    @Test
    public void asdasdas() throws Exception {
        System.out.println(isPRoot(root, (int) n, pModulus));
        for (int i = 1; i <= n; i++)
            if (powMod(root, i, pModulus) == 1)
                System.out.println("XXX");
    }

    @Test
    public void test2() throws Exception {
        long n = pow(2, 11);
        System.out.println("n: " + n);
        int c = -1;
        out:
        for (int i = 1; i < 1000; i++) {
            if (BigPrimes.isPrime(BigInteger.valueOf(i * n + 1)))
                c = i;
        }
        if (c == -1) {
            System.out.println("PIZDEC");
            return;
        }
        long prime = c * n + 1;
        System.out.println("Prime: " + prime);
        for (int i = 2; i < prime; i++) {
            if (isPRoot(i, (int) n, prime)) {
                System.out.println("root: " + i);
                break;
            }
        }
    }

    @Test
    public void test3() throws Exception {
        for (int i = 0; i <= 16; i++) {
//            System.out.println(powMod(root, i, pModulus));
            System.out.println(powMod(rootInv, i, pModulus));
        }


    }

    static long[] fft0(long[] poly, long val, long modulus) {
        MutablePolynomial p = MutablePolynomial.create(poly);
        long[] res = new long[poly.length];
        for (int i = 0; i < poly.length; i++) {
            res[i] = p.evaluate(powMod(val, i, modulus), modulus);
        }
        return res;
    }

    @Test
    public void test1() throws Exception {
        long[] poly1 = {1, 1, 0, 0, 1, 2, 0, 0, 1, 1, 1};
        long[] poly2 = {1, 0, 2, 0, 2, 0, 2, 1, 1, 0, 1};
        poly1 = Arrays.copyOf(poly1, n);
        poly2 = Arrays.copyOf(poly2, n);

        long[] fft1 = fft_3(nPower, poly1, root, pModulus);
        long[] fft2 = fft_3(nPower, poly2, root, pModulus);

        long[] fftRes = new long[fft1.length];
        for (int i = 0; i < fftRes.length; i++) {
            fftRes[i] = LongArithmetics.multiplyMod(fft1[i], fft2[i], pModulus);
        }


        long[] res = fft_3(nPower, fftRes, rootInv, pModulus);
        for (int i = 0; i < res.length; i++)
            res[i] = multiplyMod(res[i], nInv, pModulus);


        System.out.println(Arrays.toString(multiplyModClassical(poly1, 0, poly1.length, poly2, 0, poly2.length, pModulus)));
//        System.out.println(Arrays.toString(fftRes));
        System.out.println(Arrays.toString(res));
    }


    static long[] main_fft(long[] a) {
        return fft_3(nPower, a, rootPowers, pModulus);
    }

    static long[] main_fft_inv(long[] a) {
        return fft_3(nPower, a, rootInvPowers, pModulus);
    }

    static long[] fftMultiply(long[] a, long[] b) {
        int len = a.length + b.length - 1;
        a = Arrays.copyOf(a, n);
        b = Arrays.copyOf(b, n);
        long[] afft = main_fft(a);
        long[] bfft = main_fft(b);
        long[] fftRes = new long[n];

        for (int i = 0; i < fftRes.length; i++)
            fftRes[i] = LongArithmetics.multiplyMod(afft[i], bfft[i], pModulus);

        long[] res = main_fft_inv(fftRes);
        for (int i = 0; i < res.length; i++)
            res[i] = multiplyMod(res[i], nInv, pModulus);
        return Arrays.copyOf(res, len);
    }

    static long[] classicMultiply(long[] a, long[] b) {
        int aLen = a.length, bLen = b.length;
        a = Arrays.copyOf(a, n);
        b = Arrays.copyOf(b, n);

//        long[] res = MutablePolynomial.multiplyModClassical(a, 0, aLen, b, 0, bLen, pModulus);
//        return res;
        return MutablePolynomial.create(a).clone().multiply(MutablePolynomial.create(b), pModulus).data;
    }

    @Test
    public void performance_test() throws Exception {
        RandomGenerator rnd = new Well1024a();
        DescriptiveStatistics c = new DescriptiveStatistics(), f = new DescriptiveStatistics();
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                c.clear();
                f.clear();
            }
            long[] a = RandomPolynomials.randomMonicPoly(n / 2 - 100, pModulus, rnd).data;
            long[] b = RandomPolynomials.randomMonicPoly(n / 2 - 100, pModulus, rnd).data;

            long start = System.nanoTime();
            long[] classic = classicMultiply(a, b);
            long classicTime = System.nanoTime() - start;
            c.addValue(classicTime);

            start = System.nanoTime();
            long[] fft = fftMultiply(a, b);
            long fftTime = System.nanoTime() - start;
            f.addValue(fftTime);
            if (!Arrays.equals(classic, fft)) {
                System.out.println(Arrays.toString(a));
                System.out.println(Arrays.toString(b));
            }
            Assert.assertEquals(classic.length, fft.length);
            assertArrayEquals(classic, fft);
        }

        System.out.println("==== classic ==== ");
        System.out.println(c);
        System.out.println("==== fft ==== ");
        System.out.println(f);


    }


    @Test
    public void performance_test_2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        DescriptiveStatistics mul_cl = new DescriptiveStatistics(),
                fft1 = new DescriptiveStatistics(),
                fft2 = new DescriptiveStatistics(),
                product = new DescriptiveStatistics(),
                fft3 = new DescriptiveStatistics(),
                all_fft = new DescriptiveStatistics(),
                mul_fft = new DescriptiveStatistics();
        DescriptiveStatistics[] all_stats = {mul_cl, fft1, fft2, fft3, product, all_fft, mul_fft};
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                for (DescriptiveStatistics ss : all_stats) {
                    ss.clear();
                }
            }
            long[] a = RandomPolynomials.randomMonicPoly(n / 2 - 100, pModulus, rnd).data;
            long[] b = RandomPolynomials.randomMonicPoly(n / 2 - 100, pModulus, rnd).data;

            long start = System.nanoTime();
            long[] classic = classicMultiply(a, b);
            long classicTime = System.nanoTime() - start;
            mul_cl.addValue(classicTime);


            int len = a.length + b.length - 1;
            a = Arrays.copyOf(a, n);
            b = Arrays.copyOf(b, n);

            start = System.nanoTime();
            long[] afft = main_fft(a);
            long fft1time = System.nanoTime() - start;


            start = System.nanoTime();
            long[] bfft = main_fft(b);
            long fft2time = System.nanoTime() - start;

            start = System.nanoTime();
            long[] fftRes = new long[n];
            for (int k = 0; k < fftRes.length; k++)
                fftRes[k] = LongArithmetics.multiplyMod(afft[k], bfft[k], pModulus);
            long productTime = System.nanoTime() - start;


            start = System.nanoTime();
            long[] res = main_fft_inv(fftRes);
            for (int k = 0; k < res.length; k++)
                res[k] = multiplyMod(res[k], nInv, pModulus);
            long[] fft = Arrays.copyOf(res, len);
            long fft3time = System.nanoTime() - start;

            fft1.addValue(fft1time);
            fft2.addValue(fft2time);
            product.addValue(productTime);
            fft3.addValue(fft3time);
            all_fft.addValue(fft1time + fft2time + fft3time);
            mul_fft.addValue(fft1time + fft2time + fft3time + productTime);


//            if (!Arrays.equals(classic, fft)) {
//                System.out.println(Arrays.toString(a));
//                System.out.println(Arrays.toString(b));
//            }
            Assert.assertEquals(classic.length, fft.length);
            assertArrayEquals(classic, fft);
        }

        System.out.println("==== classic ==== ");
        System.out.println(mul_cl);
        System.out.println("==== fft1 ==== ");
        System.out.println(fft1);
        System.out.println("==== fft2 ==== ");
        System.out.println(fft2);
        System.out.println("==== fft3 ==== ");
        System.out.println(fft3);
        System.out.println("==== product ==== ");
        System.out.println(product);
        System.out.println("==== all_fft ==== ");
        System.out.println(all_fft);
        System.out.println("==== fft ==== ");
        System.out.println(mul_fft);


    }
}
