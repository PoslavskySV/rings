package cc.r2.core.number.qsi.factpor;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well512a;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static cc.r2.core.number.qsi.factpor.Primes.*;
import static cc.r2.core.number.qsi.factpor.Primes.QuadraticSieve;

/**
 * Created by poslavsky on 14/11/2016.
 */
public class FactorTest {

    @Test
    public void asdasdasdasv() throws Exception {
        System.out.println(factor(BigInteger.valueOf(12314324249L).multiply(BigInteger.valueOf(12312L))));
    }

    @Test
    public void as1231243() throws Exception {
        BigInteger n = BigInteger.valueOf(30011);
        n = n.multiply(BigInteger.valueOf(30103));
        n = n.multiply(BigInteger.valueOf(40111));
        n = n.multiply(BigInteger.valueOf(41113));
        n = n.multiply(BigInteger.valueOf(411101));
        n = n.multiply(BigInteger.valueOf(411113));
        n = n.multiply(BigInteger.valueOf(44111));
        System.out.println(n);
        System.out.println(QuadraticSieve(new BigInteger("18134605543"), 30000));
    }

    @Test
    public void name123() throws Exception {
        BigInteger r = new BigInteger("3124214347")
                .multiply(new BigInteger("3124214363")
                        .multiply(new BigInteger("10007")));

        System.out.println(QuadraticSieve(BigInteger.valueOf(17161), 10000));

    }

    @Test
    public void sdfsdfsdfsdf() throws Exception {
        ArrayList<BigInteger> factors = new ArrayList<>();
        BigInteger n = new BigInteger("4294967311");
        n = n.multiply(new BigInteger("4294967357"));
        n = n.multiply(new BigInteger("4294967371"));
        n = n.multiply(new BigInteger("42949679"));
        n = n.multiply(new BigInteger("42949691"));
        n = n.multiply(new BigInteger("2949701"));

        HardFactor(n, factors) ;
        BigInteger r = BigInteger.ONE;
        for (BigInteger factor : factors) {
            System.out.println(factor);
            r = factor.multiply(r);
        }

        System.out.println(n);
        System.out.println(r);
    }

    public static ArrayList<TestData> readBenchmarkingData() {
        File file = new File("/Users/poslavsky/Downloads/bench");
        if (!file.exists())
            return new ArrayList<>();
        try {
            ObjectInputStream obj = new ObjectInputStream(new FileInputStream(file));
            Object o = obj.readObject();
            return (ArrayList) o;
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    public static void writeBenchmarkingData(ArrayList<TestData> data) {
        File file = new File("/Users/poslavsky/Downloads/bench");
        if (file.exists())
            file.delete();
        try {
            ObjectOutputStream obj = new ObjectOutputStream(new FileOutputStream(file));
            obj.writeObject(data);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Test
    public void name() throws Exception {
        System.out.println(new BigInteger(512, new Random()).toString().length());
    }

    private static void addToBenchAndPrint(List<TestData> bench, TestData d) {
        System.out.println(d);
        bench.add(d);
    }

    @Test
    public void name1() throws Exception {
        List<TestData> bench = readBenchmarkingData();
        RandomGenerator rnd = new Well512a(System.currentTimeMillis());
        for (int i = 0; i < 100; i++) {
            BigInteger integer = BigInteger.valueOf(rnd.nextLong());
            if (integer.isZero())
                integer = integer.getOne();

            integer = integer.multiply(BigInteger.valueOf(rnd.nextLong()));
            if (integer.isZero())
                integer = integer.getOne();

            if (integer.compareTo(BigInteger.ZERO) < 0)
                integer = integer.negate();

            List<BigInteger> fs = TrialDivision(integer, new BigInteger("1024"));
            integer = fs.get(fs.size() - 1);

            System.out.println(integer);
            addToBenchAndPrint(bench, testPollardP1Bounded(integer, 100));
            addToBenchAndPrint(bench, testPollardP1Bounded(integer, 1000));
            addToBenchAndPrint(bench, testPollardP1Bounded(integer, 10000));
            addToBenchAndPrint(bench, testPollardP1Bounded(integer, 100000));

            addToBenchAndPrint(bench, testPollardRhoBounded(integer, 100));
            addToBenchAndPrint(bench, testPollardRhoBounded(integer, 1000));
            addToBenchAndPrint(bench, testPollardRhoBounded(integer, 10000));
            addToBenchAndPrint(bench, testPollardRhoBounded(integer, 100000));

            addToBenchAndPrint(bench, testPollardRhoRnd(integer, 1, rnd));
            addToBenchAndPrint(bench, testPollardRhoRnd(integer, 5, rnd));
            addToBenchAndPrint(bench, testPollardRhoRnd(integer, 10, rnd));
            addToBenchAndPrint(bench, testPollardRhoRnd(integer, 100, rnd));
            addToBenchAndPrint(bench, testPollardRhoRnd(integer, 1000, rnd));

//            addToBenchAndPrint(bench, testQuadraticSieve(integer, 100));
//            addToBenchAndPrint(bench, testQuadraticSieve(integer, 1000));
            addToBenchAndPrint(bench, testQuadraticSieve(integer, 10000));
            addToBenchAndPrint(bench, testQuadraticSieve(integer, 100000));
            addToBenchAndPrint(bench, testQuadraticSieve(integer, 1000000));
        }
    }

    private static TestData testQuadraticSieve(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = QuadraticSieve(num, (int) bound);
        long timing = System.nanoTime() - start;
        return new TestData("qsiev", num, factor, bound, timing);
    }

    private static TestData testPollardP1Bounded(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = PollardP1(num, bound);
        long timing = System.nanoTime() - start;
        return new TestData("rho_p", num, factor, bound, timing);
    }

    private static TestData testPollardRhoBounded(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = PollardRho(num, bound);
        long timing = System.nanoTime() - start;
        return new TestData("rho_b", num, factor, bound, timing);
    }

    private static TestData testPollardRhoRnd(BigInteger num, long tries, RandomGenerator rnd) {
        long start = System.nanoTime();
        BigInteger factor = PollardRho(num, (int) tries, rnd);
        long timing = System.nanoTime() - start;
        return new TestData("rho_r", num, factor, tries, timing);
    }

    private static class TestData implements Serializable {
        final String algorithm;
        final BigInteger num;
        final BigInteger factor;
        final long tries;
        final boolean success;
        final long timing;

        public TestData(String algorithm, BigInteger num, BigInteger factor, long tries, long timing) {
            this.algorithm = algorithm;
            this.num = num;
            this.factor = factor;
            this.tries = tries;
            this.success = factor != null;
            this.timing = timing;
        }

        @Override
        public String toString() {
            return String.format("%s\t%s\t%s\t%s\t%s", algorithm, timing, factor, tries, success);
        }
    }

    @Test
    public void qwew() throws Exception {
        System.out.println(TrialDivision(new BigInteger("295927").multiply(new BigInteger("10")),
                new BigInteger("100")));

    }
}