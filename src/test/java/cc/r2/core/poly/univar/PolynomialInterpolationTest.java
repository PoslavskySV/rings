package cc.r2.core.poly.univar;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.util.ArraysUtil;
import cc.r2.core.util.RandomUtil;
import cc.redberry.libdivide4j.FastDivision;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.univar.PolynomialInterpolation.interpolateLagrange;
import static cc.r2.core.poly.univar.PolynomialInterpolation.interpolateNewton;
import static cc.r2.core.util.TimeUnits.statisticsNanotime;
import static cc.redberry.libdivide4j.FastDivision.modSignedFast;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class PolynomialInterpolationTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        long
                modulus = 17,
                points[] = {0, 1, 2},
                values[] = {1, 2, 3};
        lMutablePolynomialZp expected = lMutablePolynomialZ.create(1, 1).modulus(modulus);
        assertEquals(expected, interpolateLagrange(points, values, modulus));
        assertEquals(expected, interpolateNewton(points, values, modulus));
    }

    @Test
    public void test2() throws Exception {
        long
                modulus = 173,
                points[] = {0, 1, 2, 4},
                values[] = {1, 2, 3, 55};
        lMutablePolynomialZp expected = lMutablePolynomialZ.create(1, 34, 37, 103).modulus(modulus);
        assertEquals(expected, interpolateLagrange(points, values, modulus));
        assertEquals(expected, interpolateNewton(points, values, modulus));
    }

    @Test
    public void test3_random_performance() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics lagrange = new DescriptiveStatistics(), newton = new DescriptiveStatistics();
        int nIterations = its(1000, 1000);
        for (int n = 0; n < nIterations; n++) {
            if (nIterations / 10 == n) {
                lagrange.clear(); newton.clear();
            }
            long[] points = RandomUtil.randomLongArray(rndd.nextInt(15, 25), 0, Short.MAX_VALUE, rnd);
            long[] values = RandomUtil.randomLongArray(points.length, 0, Short.MAX_VALUE, rnd);
            long modulus = getModulusRandom(rndd.nextInt(5, 30));

            FastDivision.Magic magic = FastDivision.magicSigned(modulus);
            for (int i = 0; i < points.length; ++i) {
                points[i] = modSignedFast(points[i], magic);
                values[i] = modSignedFast(values[i], magic);
            }
            points = ArraysUtil.getSortedDistinct(points);
            values = Arrays.copyOf(values, points.length);

            long start;
            start = System.nanoTime();
            lMutablePolynomialZp pLagrange = interpolateLagrange(points, values, modulus);
            lagrange.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            lMutablePolynomialZp pNewton = interpolateNewton(points, values, modulus);
            newton.addValue(System.nanoTime() - start);
            assertEquals(pLagrange, pNewton);
            assertInterpolation(pLagrange, points, values);
        }

        System.out.println("Lagrange : " + statisticsNanotime(lagrange));
        System.out.println("Newton   : " + statisticsNanotime(newton));
    }

    private static void assertInterpolation(lMutablePolynomialZp poly, long[] points, long[] values) {
        for (int i = 0; i < points.length; i++)
            assertEquals(values[i], poly.evaluate(points[i]));
    }
}