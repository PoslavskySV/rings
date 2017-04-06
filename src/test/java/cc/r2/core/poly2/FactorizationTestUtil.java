package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.List;

import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreePart;
import static cc.r2.core.poly2.SquareFreeFactorization.isSquareFree;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
public final class FactorizationTestUtil {
    public FactorizationTestUtil() {}

    static <T extends IMutablePolynomial<T>> void assertFactorization(T poly, FactorDecomposition<T> factorization) {
        assertEquals(poly, factorization.toPolynomial(poly));
    }

    static <T extends lMutablePolynomialAbstract<T>> void assertFactorization(T poly, long factor, List<T> factorization) {
        assertEquals(poly, factorization.stream().reduce(poly.createConstant(factor), (a, b) -> a.clone().multiply(b)));
    }

    static <T extends bMutablePolynomialAbstract<T>> void assertFactorization(T poly, BigInteger factor, List<T> factorization) {
        assertEquals(poly, factorization.stream().reduce(poly.createConstant(factor), (a, b) -> a.clone().multiply(b)));
    }

    static <T extends IMutablePolynomialZp<T>> void assertDistinctDegreeFactorization(T poly, FactorDecomposition<T> factorization) {
        for (int i = 0; i < factorization.factors.size(); i++)
            assertEquals("Factor's degree is not divisible by d.d.f. exponent",
                    0, factorization.factors.get(i).degree() % factorization.exponents.get(i));
        assertEquals(poly, factorization.toPolynomialIgnoringExponents(poly));
    }

    interface PolynomialSource {
        lMutablePolynomialZp take(long modulus);
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
        public lMutablePolynomialZp take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);

            lMutablePolynomialZp algebra = lMutablePolynomialZp.zero(modulus);
            long[] data = new long[degree + 1];
            data[0] = 1;
            for (int i = 1; i <= degree; i++)
                data[i] = algebra.addMod(algebra.multiplyMod(data[i - 1], data[i - 1]), 1);
            ArraysUtil.reverse(data, 0, data.length);

            return lMutablePolynomialZ.create(data).modulus(modulus, false);
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
        public lMutablePolynomialZp take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);
            return lMutablePolynomialZp.createMonomial(modulus, 1, degree).addMonomial(1, 1).addMonomial(1, 0);
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
        public lMutablePolynomialZp take(long modulus) {
            lMutablePolynomialZp poly;
            poly = lMutablePolynomialZ.create(rndd.nextLong(1, modulus)).modulus(modulus, false);
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

    private static abstract class AbstarctRandomSource implements PolynomialSource {
        final RandomGenerator rnd;
        final int minDegree, maxDegree;
        final boolean ensureSquareFree;
        final RandomDataGenerator rndd;

        public AbstarctRandomSource(RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            this.rnd = rnd;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.ensureSquareFree = ensureSquareFree;
            this.rndd = new RandomDataGenerator(rnd);
        }
    }

    static final class RandomSource extends AbstarctRandomSource {
        public RandomSource(RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            super(rnd, minDegree, maxDegree, ensureSquareFree);
        }

        @Override
        public lMutablePolynomialZp take(long modulus) {
            lMutablePolynomialZp poly = RandomPolynomials.randomMonicPoly(minDegree + rnd.nextInt(maxDegree - minDegree + 1), modulus, rnd).multiply(rndd.nextLong(1, modulus - 1));

            if (ensureSquareFree) {
                poly = SquareFreePart(poly);
                assert isSquareFree(poly);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }

    static final class RandomFactorableSource implements PolynomialSource {
        final int nFactors;
        final PolynomialSource pSource;
        final boolean ensureSquareFree;

        public RandomFactorableSource(int nFactors, RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            this.nFactors = nFactors;
            this.ensureSquareFree = ensureSquareFree;
            this.pSource = new RandomSource(rnd, minDegree, maxDegree, ensureSquareFree);
        }

        @Override
        public lMutablePolynomialZp take(long modulus) {
            lMutablePolynomialZp poly = lMutablePolynomialZp.one(modulus);
            for (int i = 0; i < nFactors; i++)
                poly = poly.multiply(pSource.take(modulus));

            if (ensureSquareFree) {
                poly = SquareFreePart(poly);
                assert isSquareFree(poly);
            }
            return poly;
        }
    }
}
