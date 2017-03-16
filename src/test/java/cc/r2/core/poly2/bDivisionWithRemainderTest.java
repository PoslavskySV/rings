package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly2.DivisionWithRemainder.divideAndRemainderFast;
import static cc.r2.core.poly2.bDivisionWithRemainder.divideAndRemainderClassic;
import static cc.r2.core.poly2.bDivisionWithRemainder.divideAndRemainderFast;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class bDivisionWithRemainderTest extends AbstractPolynomialTest {

    @Test
    public void test1() throws Exception {
        bMutablePolynomialMod bDividend = bMutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1).modulus(BigInteger.valueOf(7));
        MutablePolynomialMod lDividend = bDividend.toLong();


        bMutablePolynomialMod bDivider = bMutablePolynomialZ.create(1, 2, 3, 3, 2, 1).modulus(BigInteger.valueOf(7));
        MutablePolynomialMod lDivider = bDivider.toLong();

        bMutablePolynomialMod[] bqd = divideAndRemainderFast(bDividend, bDivider, true);
        MutablePolynomialMod[] lqd = divideAndRemainderFast(lDividend, lDivider, true);
        Assert.assertArrayEquals(new MutablePolynomialMod[]{bqd[0].toLong(), bqd[1].toLong()}, lqd);
    }


    @Test
    @Benchmark(runAnyway = true)
    public void test19_FastDivisionPerformance() throws Exception {
        BigInteger modulus = new BigInteger("1247842098624308285367648396697");//BigPrimes.nextPrime(new BigInteger(100, rnd));
        RandomGenerator rnd = getRandom();
        bMutablePolynomialMod divider = RandomPolynomials.randomMonicPoly(128, modulus, rnd);

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        bDivisionWithRemainder.InverseModMonomial invRev = bDivisionWithRemainder.fastDivisionPreConditioning(divider);
        long nIterations = its(1000, 5000);
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                classic.clear();
                fast.clear();
            }
            bMutablePolynomialZ dividendZ = RandomPolynomials.randomPoly(3 * divider.degree / 2, modulus, rnd);
            bMutablePolynomialMod dividend = dividendZ.modulus(modulus);

            long start = System.nanoTime();
            bMutablePolynomialMod[] qdPlain = divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            bMutablePolynomialMod[] qdNewton = divideAndRemainderFast(dividend, divider, invRev, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals(qdPlain, qdNewton);
        }

        System.out.println("==== Plain ====");
        System.out.println(classic.getPercentile(50));

        System.out.println("==== Fast ====");
        System.out.println(fast.getPercentile(50));
    }


    @Test
    @Benchmark(runAnyway = true)
    public void test20_FastDivisionPerformance() throws Exception {
        RandomGenerator rnd = getRandom();
        BigInteger modulus = new BigInteger("1247842098624308285367648396697");//BigPrimes.nextPrime(new BigInteger(100, rnd));

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        long nIterations = its(1000, 5000);
        int dividerDegree = 126;
        int dividendDegree = 256;
        for (int i = 0; i < nIterations; i++) {
            if (i * 10 == nIterations) {
                classic.clear();
                fast.clear();
            }

            bMutablePolynomialMod divider = RandomPolynomials.randomMonicPoly(dividerDegree, modulus, rnd);
            bMutablePolynomialZ dividendZ = RandomPolynomials.randomPoly(dividendDegree, modulus, rnd);
            bMutablePolynomialMod dividend = dividendZ.modulus(modulus);
            divider.multiply(BigInteger.THREE);

            long start = System.nanoTime();
            bMutablePolynomialMod[] qdPlain = divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            bMutablePolynomialMod[] qdNewton = divideAndRemainderFast(dividend, divider, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals("dividend = " + dividend + ";\ndivider = " + divider + ";\n", qdPlain, qdNewton);
        }


        System.out.println("==== Plain ====");
        System.out.println(classic.getMean());

        System.out.println("==== Fast ====");
        System.out.println(fast.getMean());
    }

    @Test
    public void test3() throws Exception {
        BigInteger modulus = new BigInteger("998427238390739620139");
        bMutablePolynomialMod dividend = bMutablePolynomialZ.parse("989441076315244786644+174683251098354358x^1+2939699558711223765x^2+993164729241539182424x^3+8652504087827847685x^4+2978039521215483585x^5+5687372540827878771x^6+3684693598277313443x^7+3034113231916032517x^8+1842720927561159970x^9+1401489172494884190x^10").modulus(modulus);
        bMutablePolynomialMod divider = bMutablePolynomialZ.parse("718119058879299323824+59748620370951943044x^1+27715597040703811206x^2+3x^3").modulus(modulus);

        bMutablePolynomialMod[] classic = divideAndRemainderClassic(dividend, divider, true);
        bMutablePolynomialMod[] fast = divideAndRemainderFast(dividend, divider, true);

        System.out.println(Arrays.toString(classic));
        System.out.println(Arrays.toString(fast));
        assertQuotientRemainder(dividend, divider, classic);
    }

    private static <T extends bMutablePolynomialAbstract<T>> void assertQuotientRemainder(T dividend, T divider, T[] qr) {
        if (qr == null) return;
        assertEquals(dividend, divider.clone().multiply(qr[0]).add(qr[1]));
    }
}