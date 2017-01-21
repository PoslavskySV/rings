package cc.r2.core.polynomial;

import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import static cc.r2.core.polynomial.RandomPolynomials.randomMonicPoly;
import static cc.r2.core.polynomial.SquareFreeFactorization.SquareFreePart;
import static cc.r2.core.polynomial.SquareFreeFactorization.isSquareFree;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
public final class FactorizationTestUtil {
    public FactorizationTestUtil() {}

    static void assertFactorization(MutablePolynomial poly, FactorDecomposition factorization) {
        assertEquals(poly, factorization.toPolynomial());
    }

    static void assertFactorization(MutablePolynomial poly, FactorDecomposition factorization, long modulus) {
        assertEquals(poly.clone().modulus(modulus), factorization.toPolynomial(modulus));
    }

    static void assertDistinctDegreeFactorization(MutablePolynomial poly, FactorDecomposition factorization, long modulus) {
        assertEquals(poly.clone().modulus(modulus), factorization.toPolynomialIgnoringExponents(modulus));
    }

    interface PolynomialSource {
        MutablePolynomial take(long modulus);
    }

    static final class WithTiming<T> {
        final T val;
        final long nanoSeconds;

        WithTiming(T val, long nanoSeconds) {
            this.val = val;
            this.nanoSeconds = nanoSeconds;
        }

        @Override
        public String toString() {
            return String.valueOf(nanoSeconds);
        }
    }

    static final class ShoupSource implements PolynomialSource {
        final int minDegree, maxDegree;
        final RandomGenerator rnd;

        public ShoupSource(RandomGenerator rnd, int minDegree, int maxDegree) {
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.rnd = rnd;
        }

        @Override
        public MutablePolynomial take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);

            long[] data = new long[degree + 1];
            data[0] = 1;
            for (int i = 1; i <= degree; i++)
                data[i] = LongArithmetics.addMod(LongArithmetics.multiplyMod(data[i - 1], data[i - 1], modulus), 1, modulus);
            ArraysUtil.reverse(data, 0, data.length);

            return MutablePolynomial.create(data);
        }
    }

    static final class GathenSource implements PolynomialSource {
        final int minDegree, maxDegree;
        final RandomGenerator rnd;

        public GathenSource(RandomGenerator rnd, int minDegree, int maxDegree) {
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.rnd = rnd;
        }

        @Override
        public MutablePolynomial take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);
            return MutablePolynomial.createMonomial(1, degree).addMonomial(1, 1).addMonomial(1, 0);
        }
    }

    static final class FactorableSource implements PolynomialSource {
        final RandomGenerator rnd;
        final int minNBase, maxNBase;
        final boolean ensureSquareFree;
        final RandomDataGenerator rndd;

        public FactorableSource(RandomGenerator rnd, int minNBase, int maxNBase, boolean ensureSquareFree) {
            this.rnd = rnd;
            this.minNBase = minNBase;
            this.maxNBase = maxNBase;
            this.ensureSquareFree = ensureSquareFree;
            this.rndd = new RandomDataGenerator(rnd);
        }

        @Override
        public MutablePolynomial take(long modulus) {
            MutablePolynomial poly;
            poly = MutablePolynomial.create(rndd.nextLong(1, modulus)).modulus(modulus);
            int nBases = rndd.nextInt(minNBase, maxNBase);
            for (int j = 1; j <= nBases; ++j)
                poly = poly.multiply(randomMonicPoly(j, modulus, rnd), modulus);

            if (ensureSquareFree) {
                poly = SquareFreePart(poly, modulus);
                assert isSquareFree(poly, modulus);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }

    static final class RandomSource implements PolynomialSource {
        final RandomGenerator rnd;
        final int minDegree, maxDegree;
        final boolean ensureSquareFree;
        final RandomDataGenerator rndd;

        public RandomSource(RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            this.rnd = rnd;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.ensureSquareFree = ensureSquareFree;
            this.rndd = new RandomDataGenerator(rnd);
        }

        @Override
        public MutablePolynomial take(long modulus) {
            MutablePolynomial poly = randomMonicPoly(minDegree + rnd.nextInt(maxDegree - minDegree + 1), modulus, rnd).multiply(rndd.nextLong(1, modulus - 1), modulus);

            if (ensureSquareFree) {
                poly = SquareFreePart(poly, modulus);
                assert isSquareFree(poly, modulus);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }
}
