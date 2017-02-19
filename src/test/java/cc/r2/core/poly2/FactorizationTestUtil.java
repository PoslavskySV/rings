package cc.r2.core.poly2;

import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreePart;
import static cc.r2.core.poly2.SquareFreeFactorization.isSquareFree;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
public final class FactorizationTestUtil {
    public FactorizationTestUtil() {}

    static <T extends MutablePolynomialAbstract<T>> void assertFactorization(T poly, FactorDecomposition<T> factorization) {
        assertEquals(poly, factorization.toPolynomial(poly));
    }

    static void assertDistinctDegreeFactorization(MutablePolynomialMod poly, FactorDecomposition<MutablePolynomialMod> factorization) {
        for (int i = 0; i < factorization.factors.size(); i++)
            assertEquals("Factor's degree is not divisible by d.d.f. exponent",
                    0, factorization.factors.get(i).degree % factorization.exponents.get(i));
        assertEquals(poly, factorization.toPolynomialIgnoringExponents(poly));
    }

    interface PolynomialSource {
        MutablePolynomialMod take(long modulus);
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
        public MutablePolynomialMod take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);

            MutablePolynomialMod algebra = MutablePolynomialMod.zero(modulus);
            long[] data = new long[degree + 1];
            data[0] = 1;
            for (int i = 1; i <= degree; i++)
                data[i] = algebra.addMod(algebra.mulMod(data[i - 1], data[i - 1]), 1);
            ArraysUtil.reverse(data, 0, data.length);

            return MutablePolynomialMod.create(modulus, data);
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
        public MutablePolynomialMod take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);
            return MutablePolynomialMod.createMonomial(modulus, 1, degree).addMonomial(1, 1).addMonomial(1, 0);
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
        public MutablePolynomialMod take(long modulus) {
            MutablePolynomialMod poly;
            poly = MutablePolynomialMod.create(modulus, rndd.nextLong(1, modulus));
            int nBases = rndd.nextInt(minNBase, maxNBase);
            for (int j = 1; j <= nBases; ++j)
                poly = poly.multiply(RandomPolynomials.randomMonicPoly(j, modulus, rnd));

            if (ensureSquareFree) {
                poly = SquareFreePart(poly);
                assert isSquareFree(poly);
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
        public MutablePolynomialMod take(long modulus) {
            MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(minDegree + rnd.nextInt(maxDegree - minDegree + 1), modulus, rnd).multiply(rndd.nextLong(1, modulus - 1));

            if (ensureSquareFree) {
                poly = SquareFreePart(poly);
                assert isSquareFree(poly);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }
}
