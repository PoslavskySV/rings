package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.IntegersZp;
import cc.r2.core.poly.MachineArithmetic;
import cc.r2.core.test.AbstractTest;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.ArraysUtil;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.Collectors;

import static cc.r2.core.poly.multivar.LinearAlgebra.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class LinearAlgebraTest extends AbstractTest {
    @Test
    public void test1() throws Exception {
        long[][] system = {{1, 2, 13}, {2, 14, 3}, {11, 2, 13}};
        long[] rhs = {1, 2, 13};
        BigInteger[] r = solve(new IntegersZp(SmallPrimes.nextPrime(12324)), convert(system), convert(rhs));
        Assert.assertArrayEquals(convert(new long[]{2467, 9625, 6865}), r);
    }

    @Test
    public void test2() throws Exception {
        long[][] lhs0 = {
                {14945154, 0, 0, 0, 0, 0, 0, 0},
                {15840518, 0, 0, 0, 0, 0, 0, 23072526},
                {0, 1, 0, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0, 0, 18556231},
                {0, 0, 7274769, 0, 0, 0, 0, 0},
                {0, 0, 1285016, 0, 0, 0, 0, 2411651},
                {0, 0, 0, 9614891, 7274769, 0, 0, 0},
                {0, 0, 0, 4514307, 1285016, 0, 0, 17488741},
                {0, 0, 0, 0, 0, 9614891, 0, 0},
                {0, 0, 0, 0, 0, 4514307, 0, 7852752},
                {0, 0, 0, 0, 0, 0, 9614891, 0},
                {0, 0, 0, 0, 0, 0, 4514307, 22089485}};
        long[] rhs0 = {1, 0, 22655879, 0, 12777324, 0, 1128298, 0, 20152010, 0, 4506067, 0};
        long modulus = 23072527;

        BigInteger[][] lhs = convert(lhs0);
        BigInteger[] rhs = convert(rhs0);

        IntegersZp domain = new IntegersZp(modulus);
        BigInteger[] solution = solve(domain, lhs, rhs);
        long[] expected = {16402965, 22655879, 11505290, 2916536, 13894224, 7600529, 2132874, 14945154};
        Assert.assertArrayEquals(convert(expected), solution);
    }

    @Test
    public void test3() throws Exception {
        long[][] lhs0 = {
                {1, 0, 0, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 5642358},
                {0, 1, 1168354, 5331039, 0, 0, 0, 0},
                {0, 1, 798805, 1341857, 0, 0, 0, 5298367},
                {0, 0, 0, 0, 3103458, 0, 0, 0}, {0, 0, 0, 0, 1168354, 0, 0, 4274594},
                {0, 0, 0, 0, 0, 1, 1168354, 0}, {0, 0, 0, 0, 0, 1, 798805, 4627257}
        };
        long[] rhs0 = {1, 0, 740880, 0, 1693671, 0, 810986, 0};
        long modulus = 5642359;
        BigInteger[] solution = new BigInteger[rhs0.length];
        SystemInfo r = solve(new IntegersZp(modulus), convert(lhs0), convert(rhs0), solution);
        Assert.assertEquals(SystemInfo.UnderDetermined, r);
    }

    @Test
    public void test3a() throws Exception {
        long[][] lhs0 = {
                {1, 0, 0, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 5642358},
                {0, 1, 1168354, 5331039, 0, 0, 0, 0},
                {0, 1, 798805, 1341857, 0, 0, 0, 5298367},
                {0, 0, 1, 0, 3103458, 0, 0, 0}, {0, 0, 0, 0, 1168354, 0, 0, 4274594},
                {0, 0, 0, 0, 0, 1, 1168354, 0}, {0, 0, 0, 0, 0, 1, 798805, 4627257}
        };
        long[] rhs0 = {1, 0, 740880, 0, 1693671, 0, 810986, 0};
        long modulus = 5642359;

//        System.out.println(prettyMatrix(lhs0));


        BigInteger[][] lhs = convert(lhs0);
        BigInteger[] rhs = convert(rhs0);

        IntegersZp domain = new IntegersZp(modulus);
        rowEchelonForm(domain, lhs, rhs);

//        System.out.println(prettyMatrix(lhs0));
//        System.out.println(prettyMatrix(lhs));
//        System.out.println(Arrays.toString(rhs));

        BigInteger[] solution = solve(domain, lhs, rhs);
        long[] expected = {1, 561035, 0, 2317604, 6, 3367849, 8000, 1};
        Assert.assertArrayEquals(convert(expected), solution);
    }

    @Test
    public void testVandermonde1() throws Exception {
        int modulus = 23;
        IntegersZp domain = new IntegersZp(modulus);
        long[] vnd = {1, 2, 3, 4};
        long[][] lhs = {
                {1, 1, 1, 1},
                {1, 2, 2 * 2, 2 * 2 * 2},
                {1, 3, 3 * 3, (3 * 3 * 3) % modulus},
                {1, 4, (4 * 4) % modulus, (4 * 4 * 4) % modulus}
        };
        long[] rhs = {11, 2, 13, 4};
        Assert.assertArrayEquals(
                solve(domain, convert(lhs), convert(rhs)),
                solveVandermonde(domain, convert(vnd), convert(rhs)));

        transposeSquare(lhs);

        Assert.assertArrayEquals(
                solve(domain, convert(lhs), convert(rhs)),
                solveVandermondeT(domain, convert(vnd), convert(rhs)));
    }

    @Test
    public void testVandermonde2() throws Exception {
        int modulus = 23;
        IntegersZp domain = new IntegersZp(modulus);
        long[] vnd = {2, 3, 4, 5};
        long[][] lhs = {
                {1, 2, 2 * 2, 2 * 2 * 2},
                {1, 3, 3 * 3, (3 * 3 * 3) % modulus},
                {1, 4, (4 * 4) % modulus, (4 * 4 * 4) % modulus},
                {1, 5, (5 * 5) % modulus, (5 * 5 * 5) % modulus}
        };
        long[] rhs = {11, 2, 13, 4};
        Assert.assertArrayEquals(
                solve(domain, convert(lhs), convert(rhs)),
                solveVandermonde(domain, convert(vnd), convert(rhs)));

        transposeSquare(lhs);

        Assert.assertArrayEquals(
                solve(domain, convert(lhs), convert(rhs)),
                solveVandermondeT(domain, convert(vnd), convert(rhs)));
    }

    @Test
    public void testVandermonde3() throws Exception {
        int modulus = 13;
        IntegersZp domain = new IntegersZp(modulus);
        long[] vnd = {2, 3, 4, 5};
        long[][] lhs = {
                {1, 2, 2 * 2, 2 * 2 * 2},
                {1, 3, 3 * 3, (3 * 3 * 3) % modulus},
                {1, 4, (4 * 4) % modulus, (4 * 4 * 4) % modulus},
                {1, 5, (5 * 5) % modulus, (5 * 5 * 5) % modulus}
        };
        long[] rhs = {11, 2, 13, 4};
        Assert.assertArrayEquals(
                solve(domain, convert(lhs), convert(rhs)),
                solveVandermonde(domain, convert(vnd), convert(rhs)));

        transposeSquare(lhs);
    }

    @Benchmark(runAnyway = true)
    @Test
    public void testVandermondePerformance() throws Exception {
        RandomGenerator rnd = getRandom();

        DescriptiveStatistics
                gauss = new DescriptiveStatistics(),
                vand = new DescriptiveStatistics(),
                gaussT = new DescriptiveStatistics(),
                vandT = new DescriptiveStatistics();
        int size = 8;
        long modulus = SmallPrimes.nextPrime(123456789);
        IntegersZp domain = new IntegersZp(modulus);
        int nIterations = 1000;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.asList(gauss, vand, gaussT, vandT).forEach(DescriptiveStatistics::clear);

            long[] k = new long[size];
            for (int i = 0; i < size; i++)
                k[i] = rnd.nextInt((int) modulus);

            long[][] lhs = new long[size][size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    lhs[i][j] = MachineArithmetic.powMod(k[i], j, modulus);


            long[] rhs = new long[size];
            for (int i = 0; i < rhs.length; i++)
                rhs[i] = rnd.nextInt((int) modulus);

            BigInteger[][] blhs = convert(lhs);
            BigInteger[] brhs = convert(rhs), vndCfs = convert(k);

            BigInteger[] gSol, vSol;
            long start;

            BigInteger[][] tmp = deepClone(blhs);
            start = System.nanoTime();
            gSol = solve(domain, tmp, brhs.clone());
            gauss.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            vSol = solveVandermonde(domain, vndCfs, brhs.clone());
            vand.addValue(System.nanoTime() - start);

            Assert.assertArrayEquals(gSol, vSol);

            transposeSquare(blhs);

            tmp = deepClone(blhs);
            start = System.nanoTime();
            gSol = solve(domain, tmp, brhs.clone());
            gaussT.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            vSol = solveVandermondeT(domain, vndCfs.clone(), brhs.clone());
            vandT.addValue(System.nanoTime() - start);

            Assert.assertArrayEquals(gSol, vSol);
        }

        System.out.println("Gauss  : " + TimeUnits.statisticsNanotime(gauss));
        System.out.println("GaussT : " + TimeUnits.statisticsNanotime(gaussT));
        System.out.println("Vandermonde  : " + TimeUnits.statisticsNanotime(vand));
        System.out.println("VandermondeT : " + TimeUnits.statisticsNanotime(vandT));
    }

    public static BigInteger[][] deepClone(BigInteger[][] input) {
        BigInteger[][] res = new BigInteger[input.length][];
        for (int i = res.length - 1; i >= 0; --i)
            res[i] = input[i].clone();
        return res;
    }

    static BigInteger[] convert(long[] arr) {
        BigInteger[] r = new BigInteger[arr.length];
        for (int i = 0; i < arr.length; i++)
            r[i] = BigInteger.valueOf(arr[i]);
        return r;
    }

    static BigInteger[][] convert(long[][] arr) {
        BigInteger[][] r = new BigInteger[arr.length][];
        for (int i = 0; i < arr.length; i++)
            r[i] = convert(arr[i]);
        return r;
    }

    private static String padding(char c, int len) {
        return new String(ArraysUtil.arrayOf(c, len));
    }

    private static String padd(String str, int newLen) {
        return padding(' ', newLen - str.length()) + str;
    }

    public static String prettyMatrix(Object[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = matrix[i][j].toString();
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        return sb.append("}").toString();
    }

    public static String prettyMatrix(long[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = Long.toString(matrix[i][j]);
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        return sb.append("}").toString();
    }
}