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

import static org.junit.Assert.assertEquals;

//import static cc.r2.core.polynomial.LongArithmetics.*;

/**
 * Created by poslavsky on 14/01/2017.
 */
@Ignore
public class FFT {

    static int nBits(long v) {
        return 64 - Long.numberOfLeadingZeros(v - 1);
    }

    static final class Magic {
        final long divisor;
        final long magic;
        final int shift;
        final int nBits;

        public Magic(long divisor, long magic, int shift) {
            this.divisor = divisor;
            this.magic = magic;
            this.shift = shift;
            this.nBits = nBits(divisor);
        }

        @Override
        public String toString() {
            return "Magic{" + "magic=" + magic + ", shift=" + shift + '}';
        }
    }

    static final Magic calculateMagics0(long d) {
        assert d < Integer.MAX_VALUE;

        long p;
        long ad, anc, delta, q1, r1, q2, r2, t;
        long two31 = 0x80000000L;     // 2**31.

        ad = Math.abs(d);
        t = two31 + (d >> 31);
        anc = t - 1 - t % ad;     // Absolute value of nc.
        p = 31;                   // Init. p.
        q1 = two31 / anc;         // Init. q1 = 2**p/|nc|.
        r1 = two31 - q1 * anc;    // Init. r1 = rem(2**p, |nc|).
        q2 = two31 / ad;          // Init. q2 = 2**p/|d|.
        r2 = two31 - q2 * ad;     // Init. r2 = rem(2**p, |d|).
        do {
            p = p + 1;
            q1 = 2 * q1;           // Update q1 = 2**p/|nc|.
            r1 = 2 * r1;           // Update r1 = rem(2**p, |nc|).
            if (r1 >= anc) {       // (Must be an unsigned
                q1 = q1 + 1;       // comparison here).
                r1 = r1 - anc;
            }
            q2 = 2 * q2;           // Update q2 = 2**p/|d|.
            r2 = 2 * r2;           // Update r2 = rem(2**p, |d|).
            if (r2 >= ad) {        // (Must be an unsigned
                q2 = q2 + 1;       // comparison here).
                r2 = r2 - ad;
            }
            delta = ad - r2;
        } while (q1 < delta || (q1 == delta && r1 == 0));

        long magM = q2 + 1;
        if (d < 0) magM = -magM; // Magic number and
        int s = (int) (p - 32);  // shift amount to return.
        //0xffffffff00000000L|
        return new Magic(d, magM, s);
    }

    static long mulhigh(long x, long y) {
        final long x_hi = x >>> 32;
        final long y_hi = y >>> 32;
        final long x_lo = x&0xFFFFFFFFL;
        final long y_lo = y&0xFFFFFFFFL;
        long result = x_lo * y_lo;
        result >>>= 32;

        result += x_hi * y_lo + x_lo * y_hi;
        result >>>= 32;
        result += x_hi * y_hi;
        return result;
    }

    static long domul(long x, long y) {
        if (x > Integer.MAX_VALUE)
            return mulhigh(x, y);
        else return x * y;
    }

    @Test
    public void asdname() throws Exception {
        System.out.println(Long.toBinaryString(0x7FFFFFFFL));

    }

    static final long remMagic(long a, Magic magic) {
        if (a - 1 > Integer.MAX_VALUE) {
            long lo = remMagic(0x7FFFFFFFL&a, magic);
            return lo + (remMagic(1L << 31, magic) * remMagic(a >>> 31, magic));
        }
        if (magic.shift == 0) {
            long q = (a * magic.magic) >>> 32;
            return a - q * magic.divisor;
        } else {
            long mag = magic.magic - (1L << 32);
            long q = (mag * a) >>> 32;
            q += a;
            q >>>= magic.shift;
            return a - q * magic.divisor;
        }
    }

    static final long remDouble(long a, long n) {
        long q = (long) ((double) a / (double) n);
        long r = a - q * n;
        if (r >= n)
            r -= n;
        else if (r < 0)
            r += n;

        return r;
    }


    static long sp_SignMask(long a) {
        return a >> (60 - 1);
    }

    static long sp_CorrectDeficit(long a, long n) {
        return a + (sp_SignMask(a)&n);
    }

    static final long remMagicDouble(long a, long n, double nInv) {
        long q = (long) (a * nInv);
        long r = a - q * n;
//        if (r >= n)
//            r -= n;
//        else if (r < 0)
//            r += n;

        r = sp_CorrectDeficit(r, n);

        return r;
    }

    @Test
    public void namedouble() throws Exception {
        long a = 1423423414213234234L;
        long n = 12312434L + Integer.MAX_VALUE;
        System.out.println(remDouble(a, n));
        System.out.println(a % n);
    }

    @Test
    public void shlag() throws Exception {
        long a = Integer.MAX_VALUE * 5L;

        System.out.println(a);
        System.out.println(0x7FFFFFFFL&a|(1L << 31) * (a >> 31));

//        System.out.println(a);
//        System.out.println(pModulus * 10000);
        //2147483647
        System.out.println(mod(a, pModulus));
        System.out.println(mod(mod(a, pModulus), pModulus));
//        System.out.println(mod(a - pModulus * 10000));
        System.out.println(a % pModulus);

    }


    @Test
    public void testmag1() throws Exception {
        for (int i = 3; i < 100000; i++) {
            if (calculateMagics0(i).shift == 0) {
                System.out.println(i);
                break;
            }
        }
        long dividend = 123143143L;
        long divisor = 33213L;
        Magic magic = calculateMagics0(divisor);
        System.out.println(magic);
        System.out.println(remMagic(dividend, magic));
        System.out.println(dividend % divisor);

    }

    static long[] reduceModBenchFast(long[] arr, Magic magic) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = remMagic(tmp[j], magic);
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    static long[] reduceModBenchPlain(long[] arr, long modulus) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = tmp[j] % modulus;
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    @Test
    public void testd7_bench() throws Exception {
//        long[] modulus = {3};
        RandomGenerator rnd = new Well1024a();
        DescriptiveStatistics plain = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        for (int i = 0; i < 100000; i++) {
            if (i == 10000) {
                fast.clear();
                plain.clear();
            }
            long[] arr = new long[1000];
            for (int j = 0; j < arr.length; j++) {
                arr[j] = rnd.nextInt();
                if (arr[j] < 0) arr[j] = -arr[j];
            }

            long modulus = rnd.nextInt();
            if (modulus < 0)
                modulus = -modulus;

            long[] f = reduceModBenchFast(arr, calculateMagics0(modulus));
            long[] p = reduceModBenchPlain(arr, modulus);

            assertEquals(f[1], p[1]);

            fast.addValue(f[0]);
            plain.addValue(p[0]);
        }

        System.out.println("==== fast ====");
        System.out.println(fast);
        System.out.println("==== plain ====");
        System.out.println(plain);
    }

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

    static long mod(long a, long b) {
        return remMagicDouble(a, b, pMagicDouble);
    }

    static long[] FastFourier0(int bits, long[] revbinData, long[] rootPowers, long modulus) {
        int n = revbinData.length;
        for (int s = 1; s <= bits; ++s) {
            int m = 1 << s;
            for (int k = 0; k <= n - 1; k += m) {
                for (int j = 0; j <= m / 2 - 1; ++j) {
                    long t = mod(rootPowers[j * n / m] * revbinData[k + j + m / 2], modulus);
                    long u = revbinData[k + j];
                    revbinData[k + j] = mod(u + t, modulus);
                    revbinData[k + j + m / 2] = mod(u - t, modulus);
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
            this.nInversed = LongArithmetics.modInverse(rootPowers.length, modulus);

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

        int nExps = 100;
        for (int i = 0; i < 1000; i++) {
            MutablePolynomial poly = RandomPolynomials.randomMonicPoly(250, pModulus, new Well1024a());
            poly = poly.modulus(pModulus);

            System.out.println(" ----- ");
            long start = System.nanoTime();
            ArrayList<MutablePolynomial> a = ExponentsNaive(poly, pModulus, nExps);
            System.out.println(System.nanoTime() - start);

            PolynomialFFT poly0 = toFFT(poly);
            InverseModMonomialFFT invMod = fastDivisionPreConditioning(poly0);

            start = System.nanoTime();
            ArrayList<PolynomialFFT> b = ExponentsNaiveFFT(poly0, invMod, pModulus, nExps);
            System.out.println(System.nanoTime() - start);

            for (int j = 0; j < nExps; j++) {
                assertEquals(a.get(j), b.get(j).toPoly());
            }
        }


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
    static final long nInv = LongArithmetics.modInverse(n, pModulus);


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


    static final Magic pMagic = calculateMagics0(pModulus);
    static final double pMagicDouble = 1.0 / pModulus;
    static final long[] rootPowers;
    static final long[] rootInvPowers;

    static {
        rootPowers = new long[n];
        rootInvPowers = new long[n];
        for (int i = 0; i < n; i++) {
            rootPowers[i] = LongArithmetics.powMod(root, i, pModulus);
            rootInvPowers[i] = LongArithmetics.powMod(rootInv, i, pModulus);
        }
    }

    static boolean isPRoot(long root, int exp, long modulus) {
        if (LongArithmetics.powMod(root, exp, modulus) != 1)
            return false;
        for (int i = 1; i < exp; i++)
            if (LongArithmetics.powMod(root, i, modulus) == 1)
                return false;
        return true;
    }

    @Test
    public void asdasdas() throws Exception {
        System.out.println(isPRoot(root, (int) n, pModulus));
        for (int i = 1; i <= n; i++)
            if (LongArithmetics.powMod(root, i, pModulus) == 1)
                System.out.println("XXX");
    }

    @Test
    public void test2() throws Exception {
        long n = LongArithmetics.pow(2, 11);
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
}
