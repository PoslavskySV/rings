package cc.r2.core.number.qsi.factpor;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well512a;

import java.io.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

import static cc.r2.core.number.qsi.factpor.Primes.*;
import static cc.r2.core.number.qsi.factpor.Primes.QuadraticSieve;

/**
 * Created by poslavsky on 14/11/2016.
 */
public class Benchmark {
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

    private static final Object mutex = new Object();

    private static volatile int warmUp = 100;

    private static void addToBenchAndPrint(TestData d) {
        synchronized (mutex) {
            if (warmUp > 0) {
                --warmUp;
                if (warmUp == 0)
                    System.out.println("JVM warmed.");
                return;
            }

            stream.println(d);
        }
    }

    public static PrintStream stream;

    public static void main(String[] args) throws Exception {
        int threads = Integer.parseInt(args[0]);
        int nTImes = Integer.parseInt(args[1]);
        stream = new PrintStream(new FileOutputStream(args[2]));

        Thread[] thread = new Thread[threads];
        for (int i = 0; i < thread.length; i++) {
            thread[i] = new Thread(new Worker(nTImes, i));
            thread[i].start();
        }

        for (int i = 0; i < thread.length; i++)
            thread[i].join();
    }

    private static class Worker implements Runnable {
        private final int nTimes;
        private final int nThread;

        public Worker(int nTimes, int nThread) {
            this.nTimes = nTimes;
            this.nThread = nThread;
        }

        @Override
        public void run() {
            RandomGenerator rnd = new Well512a(new SecureRandom().nextLong());
            for (int i = 0; i < nTimes; i++) {
                if (i % 30 == 0)
                    System.out.println(nThread + ": " + i);
                BigInteger integer = null;
                if (i % 3 == 0) {
                    //smaller ints
                    integer = new BigInteger("9999").add(BigInteger.valueOf(rnd.nextInt()));
                } else if (i % 8 == 0) {
                    //even more smaller ints
                    integer = new BigInteger("10").add(BigInteger.valueOf(rnd.nextInt(1048576)));
                } else if (i % 5 == 0) {
                    //even more smaller ints
                    integer = BigInteger.valueOf(rnd.nextInt());
                    if (integer.isZero())
                        integer = integer.getOne();

                    integer = integer.multiply(BigInteger.valueOf(rnd.nextInt()));
                    if (integer.isZero())
                        integer = integer.getOne();
                } else {
                    integer = BigInteger.valueOf(rnd.nextLong());
                    if (integer.isZero())
                        integer = integer.getOne();

                    integer = integer.multiply(BigInteger.valueOf(rnd.nextLong()));
                    if (integer.isZero())
                        integer = integer.getOne();

                }

                if (integer.isZero()) {
                    --i;
                    continue;
                }

                if (integer.compareTo(BigInteger.ZERO) < 0)
                    integer = integer.negate();

                List<BigInteger> fs = TrialDivision(integer, new BigInteger("1024"));
                integer = fs.get(fs.size() - 1);

//                addToBenchAndPrint(testFernmat(integer, 100));
//                addToBenchAndPrint(testFernmat(integer, 1000));
//                addToBenchAndPrint(testFernmat(integer, 10000));
//                addToBenchAndPrint(testFernmat(integer, 100000));

                addToBenchAndPrint(testPollardP1Bounded(integer, 100));
                addToBenchAndPrint(testPollardP1Bounded(integer, 1000));
                addToBenchAndPrint(testPollardP1Bounded(integer, 10000));
                addToBenchAndPrint(testPollardP1Bounded(integer, 100000));

                addToBenchAndPrint(testPollardRhoBounded(integer, 100));
                addToBenchAndPrint(testPollardRhoBounded(integer, 1000));
                addToBenchAndPrint(testPollardRhoBounded(integer, 10000));
                addToBenchAndPrint(testPollardRhoBounded(integer, 100000));

                addToBenchAndPrint(testPollardRhoRnd(integer, 1, rnd));
                addToBenchAndPrint(testPollardRhoRnd(integer, 5, rnd));
                addToBenchAndPrint(testPollardRhoRnd(integer, 10, rnd));
                addToBenchAndPrint(testPollardRhoRnd(integer, 100, rnd));
                addToBenchAndPrint(testPollardRhoRnd(integer, 1000, rnd));

//            addToBenchAndPrint(bench, testQuadraticSieve(integer, 100));
//            addToBenchAndPrint(bench, testQuadraticSieve(integer, 1000));
                addToBenchAndPrint(testQuadraticSieve(integer, 10000));
                addToBenchAndPrint(testQuadraticSieve(integer, 100000));
                addToBenchAndPrint(testQuadraticSieve(integer, 1000000));
            }
        }
    }

    private static TestData testFernmat(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = fermat(num, (int) bound);
        long timing = System.nanoTime() - start;
        return new TestData("fermat", num, factor, bound, timing);
    }

    private static TestData testQuadraticSieve(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = QuadraticSieve(num, (int) bound);
        long timing = System.nanoTime() - start;
        return new TestData("quadratic_sieve", num, factor, bound, timing);
    }

    private static TestData testPollardP1Bounded(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = PollardP1(num, bound);
        long timing = System.nanoTime() - start;
        return new TestData("pollard_p1", num, factor, bound, timing);
    }

    private static TestData testPollardRhoBounded(BigInteger num, long bound) {
        long start = System.nanoTime();
        BigInteger factor = PollardRho(num, bound);
        long timing = System.nanoTime() - start;
        return new TestData("pollard_rho_bounded", num, factor, bound, timing);
    }

    private static TestData testPollardRhoRnd(BigInteger num, long tries, RandomGenerator rnd) {
        long start = System.nanoTime();
        BigInteger factor = PollardRho(num, (int) tries, rnd);
        long timing = System.nanoTime() - start;
        return new TestData("pollard_rho_random", num, factor, tries, timing);
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
            return String.format("%s\t%s\t%s\t%s\t%s\t%s", algorithm, timing, num, factor, tries, success);
        }
    }

}
