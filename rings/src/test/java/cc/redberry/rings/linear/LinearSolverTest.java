package cc.redberry.rings.linear;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.Collectors;

import static cc.redberry.rings.linear.LinearSolver.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class LinearSolverTest extends AbstractTest {
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
    public void test2a() throws Exception {
        long[][] lhs0 = {
                {14945154, 0, 0, 0, 0, 0, 0, 0},
                {14945154, 0, 0, 0, 0, 0, 0, 0},
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
        long[] rhs0 = {1, 1, 1, 0, 22655879, 0, 12777324, 0, 1128298, 0, 20152010, 0, 4506067, 0};
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
        assertEquals(SystemInfo.UnderDetermined, r);
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
    public void test4() throws Exception {
        long[][] lhs = {
                {1, 2, 3, 4, 5, 6},
                {1, 2, 3, 4, 5, 6},
                {3, 4, 2, 1, 3, 2}
        };
        long[] rhs = {1, 1, 2};
        long[] solution = new long[lhs[0].length];
        IntegersZp64 ring = new IntegersZp64(17);
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test5() {
        long[][] lhs = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 16, 0, 0, 0, 0, 0, 0, 0, 2, 0, 16, 0, 1, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 1, 0, 1, 1, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 1, 0, 1, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 1, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        long[] rhs = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        long[] solution = new long[lhs[0].length];
        IntegersZp64 ring = new IntegersZp64(17);
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test6_random() {
        RandomGenerator rnd = getRandom();
        int nIterations = 1_000;
        for (long modulus : Arrays.asList(2, 17, SmallPrimes.nextPrime(1000))) {
            IntegersZp64 ring = new IntegersZp64(modulus);
            int[] nn = {10, 15, 20};
            for (boolean consistent : Arrays.asList(true, false))
                for (int nRows : nn)
                    for (int nColumns : nn)
                        testRandom(ring, nRows, nColumns, nIterations, consistent, rnd);
        }
    }

    @Test
    public void test7() {
        IntegersZp64 ring = new IntegersZp64(SmallPrimes.nextPrime(1000));
        long[][] lhs = {
                {529, 862, 507},
                {427, 116, 939},
                {281, 3, 477},
        };
        long[] rhs = {406, 346, 754};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test8() {
        IntegersZp64 ring = new IntegersZp64(SmallPrimes.nextPrime(1000));
        long[][] lhs = {
                {927, 433, 167, 418},
                {671, 760, 276, 769},
                {485, 535, 316, 204},
        };
        long[] rhs = {437, 433, 646};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test9() {
        IntegersZp64 ring = new IntegersZp64(SmallPrimes.nextPrime(1000));
        long[][] lhs = {
                {86, 300, 654},
                {421, 15, 279},
                {143, 380, 386},
        };
        long[] rhs = {701, ring.modulus(-1), 499};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test10() {
        IntegersZp64 ring = new IntegersZp64(SmallPrimes.nextPrime(1000));
        long[][] lhs = {
                {815, 0, 342},
                {170, 536, 982},
                {267, 31, 878}
        };
        long[] rhs = {433, 263, 439};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Inconsistent, r);
    }

    @Test
    public void test11() {
        IntegersZp64 ring = new IntegersZp64(17);
        long[][] lhs = {
                {5, 9, 8, 4, 0, 10, 14, 12, 16, 16},
                {2, 11, 6, 2, 9, 16, 11, 2, 4, 15},
                {15, 11, 16, 5, 9, 16, 12, 15, 2, 6},
                {0, 10, 11, 2, 0, 13, 11, 15, 16, 14},
                {8, 6, 12, 7, 2, 2, 14, 14, 6, 16},
                {16, 3, 1, 2, 3, 8, 10, 13, 3, 9},
                {1, 2, 12, 16, 1, 12, 14, 15, 10, 2},
                {6, 5, 5, 15, 11, 3, 15, 15, 3, 14},
                {1, 12, 11, 4, 1, 0, 7, 7, 11, 0},
                {14, 14, 3, 0, 15, 10, 1, 2, 10, 14}
        };
        long[] rhs = {2, 5, 13, 1, 15, 13, 12, 14, 6, 8};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Inconsistent, r);
    }

    @Test
    public void test12() {
        IntegersZp64 ring = new IntegersZp64(17);
        long[][] lhs = {
                {8, 9, 6, 13, 4, 8, 11, 12, 15, 9},
                {11, 2, 3, 5, 8, 9, 6, 2, 8, 15},
                {15, 15, 10, 6, 5, 8, 3, 13, 12, 16},
                {15, 2, 14, 14, 0, 5, 5, 0, 9, 15},
                {12, 11, 9, 9, 5, 4, 2, 2, 0, 3},
                {12, 9, 14, 5, 9, 10, 7, 4, 15, 4},
                {2, 16, 16, 5, 16, 10, 3, 0, 15, 3},
                {14, 13, 13, 0, 4, 14, 12, 13, 1, 3},
                {0, 15, 1, 10, 13, 2, 3, 4, 13, 6},
                {6, 7, 2, 11, 14, 1, 6, 6, 7, 5}
        };
        long[] rhs = {6, 10, 5, 2, 8, 8, 8, 14, 8, 5};
        long[] solution = new long[lhs[0].length];
        SystemInfo r = solve(ring, deepClone(lhs), rhs.clone(), solution, true);
        assertEquals(SystemInfo.Consistent, r);
        assertTrue(isSolution(ring, lhs, solution, rhs));
    }

    @Test
    public void test13_random() {
        RandomGenerator rnd = getRandom();
        int nIterations = 100;
        for (long modulus : Arrays.asList(2, 17, SmallPrimes.nextPrime(1000))) {
            IntegersZp ring = new IntegersZp(modulus);
            int[] nn = {10, 15, 20};
            for (boolean consistent : Arrays.asList(true, false))
                for (int nRows : nn)
                    for (int nColumns : nn)
                        testRandom(ring, nRows, nColumns, nIterations, consistent, rnd);
        }
    }


    private static void testRandom(IntegersZp64 ring,
                                   int nRows,
                                   int nColumns,
                                   int nIterations,
                                   boolean consistent,
                                   RandomGenerator rnd) {
        for (int iter = 0; iter < nIterations; ++iter) {
            long[][] lhs = new long[nRows][nColumns];
            long[] solution = new long[nColumns];

            for (int i = 0; i < Math.min(nRows, nColumns); ++i)
                solution[i] = ring.randomElement(rnd);

            for (int i = 0; i < nRows; ++i)
                for (int j = 0; j < nColumns; ++j)
                    lhs[i][j] = ring.randomElement(rnd);

            long[] rhs = multiply(ring, lhs, solution);
            assert isSolution(ring, lhs, solution, rhs);

            // rearrange rows
            for (int n = 0; n < 2 * nRows; ++n) {
                int i = rnd.nextInt(nRows), j = rnd.nextInt(nRows);
                if (i == j)
                    continue;
                long
                        f1 = ring.randomNonZeroElement(rnd),
                        f2 = ring.randomNonZeroElement(rnd);
                long[] p = add(ring, multiply(ring, lhs[i], f1), multiply(ring, lhs[j], f2));
                long[] q = add(ring, multiply(ring, lhs[i], f1), multiply(ring, lhs[j], ring.negate(f2)));

                lhs[i] = p;
                lhs[j] = q;

                long rp = ring.add(ring.multiply(rhs[i], f1), ring.multiply(rhs[j], f2));
                long rq = ring.add(ring.multiply(rhs[i], f1), ring.multiply(rhs[j], ring.negate(f2)));

                rhs[i] = rp;
                rhs[j] = rq;

                if (!consistent) {
                    rhs[j] = ring.add(rhs[j], 1);
                    rhs[i] = ring.subtract(rhs[i], 1);
                }
            }

            // rearrange columns
            for (int n = 0; n < nColumns; ++n) {
                int i = rnd.nextInt(nColumns), j = rnd.nextInt(nColumns);
                swapColumns(lhs, i, j);
                ArraysUtil.swap(solution, i, j);
            }

            try {
                long[] actualSolution = new long[solution.length];
                SystemInfo solve = LinearSolver.solve(ring, deepClone(lhs), rhs.clone(), actualSolution, true);
                if (consistent) {
                    assertTrue(solve == SystemInfo.Consistent);
                    assertTrue(isSolution(ring, lhs, actualSolution, rhs));
                } else
                    assertTrue(solve == SystemInfo.Inconsistent || isSolution(ring, lhs, actualSolution, rhs));
            } catch (Throwable t) {
                System.out.println(prettyMatrix(lhs));
                System.out.println(Arrays.toString(rhs));
                throw new RuntimeException(t);
            }
        }
    }

    private static void swapColumns(long[][] matrix, int i, int j) {
        if (i == j)
            return;
        for (long[] aMatrix : matrix)
            ArraysUtil.swap(aMatrix, i, j);
    }

    private static long[] multiply(IntegersZp64 ring, long[][] matrix, long[] column) {
        long[] result = new long[matrix.length];
        for (int i = 0; i < matrix.length; ++i)
            for (int j = 0; j < matrix[i].length; ++j)
                result[i] = ring.add(result[i], ring.multiply(matrix[i][j], column[j]));
        return result;
    }

    private static long[] add(IntegersZp64 ring, long[] a, long[] b) {
        long[] result = new long[a.length];
        for (int i = 0; i < a.length; ++i)
            result[i] = ring.add(result[i], ring.add(a[i], b[i]));
        return result;
    }

    private static long[] multiply(IntegersZp64 ring, long[] a, long b) {
        long[] result = new long[a.length];
        for (int i = 0; i < a.length; ++i)
            result[i] = ring.multiply(a[i], b);
        return result;
    }


    private static <E> void testRandom(Ring<E> ring,
                                       int nRows,
                                       int nColumns,
                                       int nIterations,
                                       boolean consistent,
                                       RandomGenerator rnd) {
        for (int iter = 0; iter < nIterations; ++iter) {
            E[][] lhs = ring.createZeroesArray2d(nRows, nColumns);
            E[] solution = ring.createZeroesArray(nColumns);

            for (int i = 0; i < Math.min(nRows, nColumns); ++i)
                solution[i] = ring.randomElement(rnd);

            for (int i = 0; i < nRows; ++i)
                for (int j = 0; j < nColumns; ++j)
                    lhs[i][j] = ring.randomElement(rnd);

            E[] rhs = multiply(ring, lhs, solution);
            assert isSolution(ring, lhs, solution, rhs);

            // rearrange rows
            for (int n = 0; n < 2 * nRows; ++n) {
                int i = rnd.nextInt(nRows), j = rnd.nextInt(nRows);
                if (i == j)
                    continue;
                E
                        f1 = ring.randomNonZeroElement(rnd),
                        f2 = ring.randomNonZeroElement(rnd);
                E[] p = add(ring, multiply(ring, lhs[i], f1), multiply(ring, lhs[j], f2));
                E[] q = add(ring, multiply(ring, lhs[i], f1), multiply(ring, lhs[j], ring.negate(f2)));

                lhs[i] = p;
                lhs[j] = q;

                E rp = ring.add(ring.multiply(rhs[i], f1), ring.multiply(rhs[j], f2));
                E rq = ring.add(ring.multiply(rhs[i], f1), ring.multiply(rhs[j], ring.negate(f2)));

                rhs[i] = rp;
                rhs[j] = rq;

                if (!consistent) {
                    rhs[j] = ring.add(rhs[j], ring.getOne());
                    rhs[i] = ring.subtract(rhs[i], ring.getOne());
                }
            }

            // rearrange columns
            for (int n = 0; n < nColumns; ++n) {
                int i = rnd.nextInt(nColumns), j = rnd.nextInt(nColumns);
                swapColumns(lhs, i, j);
                ArraysUtil.swap(solution, i, j);
            }

            try {
                E[] actualSolution = ring.createArray(solution.length);
                SystemInfo solve = LinearSolver.solve(ring, deepClone(lhs), rhs.clone(), actualSolution, true);
                if (consistent) {
                    assertTrue(solve == SystemInfo.Consistent);
                    assertTrue(isSolution(ring, lhs, actualSolution, rhs));
                } else
                    assertTrue(solve == SystemInfo.Inconsistent || isSolution(ring, lhs, actualSolution, rhs));
            } catch (Throwable t) {
                System.out.println(prettyMatrix(lhs));
                System.out.println(Arrays.toString(rhs));
                throw new RuntimeException(t);
            }
        }
    }

    private static void swapColumns(Object[][] matrix, int i, int j) {
        if (i == j)
            return;
        for (Object[] aMatrix : matrix)
            ArraysUtil.swap(aMatrix, i, j);
    }

    private static <E> E[] multiply(Ring<E> ring, E[][] matrix, E[] column) {
        E[] result = ring.createZeroesArray(matrix.length);
        for (int i = 0; i < matrix.length; ++i)
            for (int j = 0; j < matrix[i].length; ++j)
                result[i] = ring.add(result[i], ring.multiply(matrix[i][j], column[j]));
        return result;
    }

    private static <E> E[] add(Ring<E> ring, E[] a, E[] b) {
        E[] result = ring.createZeroesArray(a.length);
        for (int i = 0; i < a.length; ++i)
            result[i] = ring.add(result[i], ring.add(a[i], b[i]));
        return result;
    }

    private static <E> E[] multiply(Ring<E> ring, E[] a, E b) {
        E[] result = ring.createArray(a.length);
        for (int i = 0; i < a.length; ++i)
            result[i] = ring.multiply(a[i], b);
        return result;
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

    static boolean isSolution(IntegersZp64 ring, long[][] lhs, long[] solution, long[] rhs) {
        return Arrays.equals(multiply(ring, lhs, solution), rhs);
    }

    static <E> boolean isSolution(Ring<E> ring, E[][] lhs, E[] solution, E[] rhs) {
        return Arrays.equals(multiply(ring, lhs, solution), rhs);
    }

    static BigInteger[][] deepClone(BigInteger[][] input) {
        BigInteger[][] res = new BigInteger[input.length][];
        for (int i = res.length - 1; i >= 0; --i)
            res[i] = input[i].clone();
        return res;
    }

    static <E> E[][] deepClone(E[][] input) {
        E[][] res = input.clone();
        for (int i = res.length - 1; i >= 0; --i)
            res[i] = input[i].clone();
        return res;
    }

    static long[][] deepClone(long[][] input) {
        long[][] res = new long[input.length][];
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