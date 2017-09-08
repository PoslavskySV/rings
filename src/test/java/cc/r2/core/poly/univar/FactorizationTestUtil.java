package cc.r2.core.poly.univar;

import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
final class FactorizationTestUtil {
    public FactorizationTestUtil() {}

    public static <T extends IUnivariatePolynomial<T>> void assertDistinctDegreeFactorization(T poly, FactorDecomposition<T> factorization) {
        for (int i = 0; i < factorization.factors.size(); i++)
            assertEquals("Factor's degree is not divisible by d.d.f. exponent",
                    0, factorization.factors.get(i).degree() % factorization.exponents.get(i));
        assertEquals(poly, factorization.toPolynomialIgnoringExponents());
    }


    public static <T extends lUnivariatePolynomialAbstract<T>> void assertFactorization(T poly, long factor, List<T> factorization) {
        assertEquals(poly, factorization.stream().reduce(poly.createConstant(factor), (a, b) -> a.clone().multiply(b)));
    }

    public static <E> void assertFactorization(UnivariatePolynomial<E> poly, E factor, List<UnivariatePolynomial<E>> factorization) {
        assertEquals(poly, factorization.stream().reduce(poly.createConstant(factor), (a, b) -> a.clone().multiply(b)));
    }

    public interface PolynomialSource {
        lUnivariatePolynomialZp take(long modulus);
    }

    public static final class WithTiming<T> {
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

    public static final class ShoupSource implements PolynomialSource {
        final int minDegree, maxDegree;
        final RandomGenerator rnd;

        public ShoupSource(RandomGenerator rnd, int minDegree, int maxDegree) {
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.rnd = rnd;
        }

        @Override
        public lUnivariatePolynomialZp take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);

            lUnivariatePolynomialZp algebra = lUnivariatePolynomialZp.zero(modulus);
            long[] data = new long[degree + 1];
            data[0] = 1;
            for (int i = 1; i <= degree; i++)
                data[i] = algebra.add(algebra.multiply(data[i - 1], data[i - 1]), 1);
            ArraysUtil.reverse(data, 0, data.length);

            return lUnivariatePolynomialZ.create(data).modulus(modulus, false);
        }
    }

    public static final class GathenSource implements PolynomialSource {
        final int minDegree, maxDegree;
        final RandomGenerator rnd;

        public GathenSource(RandomGenerator rnd, int minDegree, int maxDegree) {
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.rnd = rnd;
        }

        @Override
        public lUnivariatePolynomialZp take(long modulus) {
            int degree = minDegree + rnd.nextInt(maxDegree - minDegree + 1);
            return lUnivariatePolynomialZp.monomial(modulus, 1, degree).addMonomial(1, 1).addMonomial(1, 0);
        }
    }

    public static final class FactorableSource implements PolynomialSource {
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
        public lUnivariatePolynomialZp take(long modulus) {
            lUnivariatePolynomialZp poly;
            poly = lUnivariatePolynomialZ.create(rndd.nextLong(1, modulus)).modulus(modulus, false);
            int nBases = rndd.nextInt(minNBase, maxNBase);
            for (int j = 1; j <= nBases; ++j)
                poly = poly.multiply(RandomUnivariatePolynomials.randomMonicPoly(j, modulus, rnd));

            if (ensureSquareFree) {
                poly = UnivariateSquareFreeFactorization.SquareFreePart(poly);
                assert UnivariateSquareFreeFactorization.isSquareFree(poly);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }

    private static abstract class AbstractRandomSource implements PolynomialSource {
        final RandomGenerator rnd;
        final int minDegree, maxDegree;
        final boolean ensureSquareFree;
        final RandomDataGenerator rndd;

        public AbstractRandomSource(RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            this.rnd = rnd;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.ensureSquareFree = ensureSquareFree;
            this.rndd = new RandomDataGenerator(rnd);
        }
    }

    public static final class RandomSource extends AbstractRandomSource {
        public RandomSource(RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            super(rnd, minDegree, maxDegree, ensureSquareFree);
        }

        @Override
        public lUnivariatePolynomialZp take(long modulus) {
            lUnivariatePolynomialZp poly = RandomUnivariatePolynomials.randomMonicPoly(minDegree + rnd.nextInt(maxDegree - minDegree + 1), modulus, rnd).multiply(rndd.nextLong(1, modulus - 1));

            if (ensureSquareFree) {
                poly = UnivariateSquareFreeFactorization.SquareFreePart(poly);
                assert UnivariateSquareFreeFactorization.isSquareFree(poly);
            }

            if (poly.isConstant())
                return take(modulus);

            return poly;
        }
    }

    public static final class RandomFactorableSource implements PolynomialSource {
        final int nFactors;
        final PolynomialSource pSource;
        final boolean ensureSquareFree;

        public RandomFactorableSource(int nFactors, RandomGenerator rnd, int minDegree, int maxDegree, boolean ensureSquareFree) {
            this.nFactors = nFactors;
            this.ensureSquareFree = ensureSquareFree;
            this.pSource = new RandomSource(rnd, minDegree, maxDegree, ensureSquareFree);
        }

        @Override
        public lUnivariatePolynomialZp take(long modulus) {
            lUnivariatePolynomialZp poly = lUnivariatePolynomialZp.one(modulus);
            for (int i = 0; i < nFactors; i++)
                poly = poly.multiply(pSource.take(modulus));

            if (ensureSquareFree) {
                poly = UnivariateSquareFreeFactorization.SquareFreePart(poly);
                assert UnivariateSquareFreeFactorization.isSquareFree(poly);
            }
            return poly;
        }
    }
}
