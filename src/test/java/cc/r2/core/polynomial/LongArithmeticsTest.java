package cc.r2.core.polynomial;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.PrimesIterator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import static cc.r2.core.polynomial.LongArithmetics.*;
import static java.lang.Math.floorMod;
import static org.junit.Assert.*;

public class LongArithmeticsTest {

    @Test
    public void test1() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 1000; i++) {
            long a = 1 + rnd.nextInt(1000), b = 1 + rnd.nextInt(100000);
            long[] gcd = gcdExtended(a, b);
            assertEquals(0, a % gcd[0]);
            assertEquals(0, b % gcd[0]);
            assertEquals(gcd[0], gcd[1] * a + gcd[2] * b);
        }
    }

    @Test
    public void test2() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 1000; i++) {
            long aa = 1 + rnd.nextInt(1000), b = 1 + rnd.nextInt(100000);
            if (gcd(aa, b) != 1) {
                --i;
                continue;
            }
            for (long a : new long[]{aa, -aa}) {
                long l = modInverse(a, b);
                assertEquals(1L, floorMod(a * l, b));
                assertEquals(BigInteger.valueOf(a).modInverse(BigInteger.valueOf(b)).longValue(), l);
            }
        }
    }

    @Test
    public void test3() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            int x = rnd.nextInt();
            int y = rnd.nextInt();
            assertTrue(safeMultiply(x, y));
            assertTrue(safeMultiply(x, -y));
            assertTrue(safeMultiply(-x, y));
            assertTrue(safeMultiply(-x, -y));
        }

        for (int i = 0; i < 100; i++) {
            long x = ((long) Integer.MAX_VALUE) * 2 + 10 + rnd.nextInt(100);
            long y = ((long) Integer.MAX_VALUE) * 2 + 10 + rnd.nextInt(100);
            assertFalse(safeMultiply(x, y));
            assertFalse(safeMultiply(x, -y));
            assertFalse(safeMultiply(-x, y));
            assertFalse(safeMultiply(-x, -y));
        }
    }

    @Test
    public void test4() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            long x = rnd.nextLong() / 2;
            long y = rnd.nextLong() / 2;
            assertTrue(safeAdd(x, y));
            assertTrue(safeAdd(x, -y));
            assertTrue(safeAdd(-x, y));
            assertTrue(safeAdd(-x, -y));
        }

        for (int i = 0; i < 100; i++) {
            long x = Long.MAX_VALUE - 1000 + rnd.nextInt(1000);
            long y = Long.MAX_VALUE - 1000 + rnd.nextInt(1000);
            assertFalse(safeAdd(x, y));
            assertFalse(safeAdd(-x, -y));
            assertTrue(safeAdd(x, -y));
            assertTrue(safeAdd(-x, y));
        }
    }

    @Test
    public void test5() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            long x = rnd.nextLong() / 2;
            long y = rnd.nextLong() / 2;
            assertTrue(safeSubtract(x, y));
            assertTrue(safeSubtract(x, -y));
            assertTrue(safeSubtract(-x, y));
            assertTrue(safeSubtract(-x, -y));
        }

        for (int i = 0; i < 100; i++) {
            long x = Long.MAX_VALUE - 1000 + rnd.nextInt(1000);
            long y = Long.MAX_VALUE - 1000 + rnd.nextInt(1000);
            assertFalse(safeSubtract(x, -y));
            assertFalse(safeSubtract(-x, y));
            assertTrue(safeSubtract(x, y));
            assertTrue(safeSubtract(-x, -y));
        }
    }

    @Test
    public void test6() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 1000; i++) {
            int x = rnd.nextInt(1000);
            int y = rnd.nextInt(6);
//            assertEquals(pow(x, y), pow(x, y));
//            assertEquals(pow(-x, y), pow(-x, y));
            assertTrue(safePow(x, y));
            assertTrue(safePow(-x, y));
        }

        for (int i = 0; i < 1000; i++) {
            int x = 1000 + rnd.nextInt(1000);
            int y = 8 + rnd.nextInt(6);
            assertFalse(safePow(x, y));
            assertFalse(safePow(-x, y));
            try { pow(x, y);} catch (Exception e) {assertEquals(ArithmeticException.class, e.getClass());}
            try { pow(-x, y);} catch (Exception e) {assertEquals(ArithmeticException.class, e.getClass());}
        }
    }

    @Test
    public void test7() throws Exception {
        Well1024a rnd = new Well1024a();
        for (int i = 0; i < 1000; i++) {
            long a = 1 + rnd.nextInt(1000), b = 1 + rnd.nextInt(100000);
            if (gcd(a, b) != 1) {
                --i;
                continue;
            }
            long e = 1 + rnd.nextInt(100);
            assertEquals(BigInteger.valueOf(a).modPow(BigInteger.valueOf(e), BigInteger.valueOf(b)).longValue(), powMod(a, e, b));
            e = -e;
            assertEquals(BigInteger.valueOf(a).modPow(BigInteger.valueOf(e), BigInteger.valueOf(b)).longValue(), powMod(a, e, b));
        }
    }

    @Test
    public void test8() throws Exception {
        assertEquals(BigInteger.valueOf(-1).modInverse(BigInteger.valueOf(11)).longValue(), modInverse(-1, 11));
    }

    @Test
    public void test9() throws Exception {
        assertEquals(0, symMod(0, 7));
        assertEquals(3, symMod(3, 7));
        assertEquals(-3, symMod(4, 7));
    }

    @Test
    public void test10() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            long prime = new PrimesIterator(rnd.nextInt(10000)).take();
            for (int j = 0; j < 10; j++) {
                long k = rnd.nextInt();
                long symMod = symMod(k, prime);
                assertEquals(floorMod(symMod, prime), floorMod(k, prime));
            }
        }
    }
}