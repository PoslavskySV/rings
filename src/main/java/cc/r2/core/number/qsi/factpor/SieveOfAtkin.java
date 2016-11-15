package cc.r2.core.number.qsi.factpor;

import org.apache.commons.math3.primes.*;

import java.util.Arrays;
import java.util.BitSet;

/**
 * Created by poslavsky on 14/11/2016.
 */
public class SieveOfAtkin {
    private static int limit = 100000;
    private static BitSet sieve = new BitSet(limit + 1);
    private static int limitSqrt = (int) Math.sqrt((double) limit);

    public static void main(String[] args) {
        // there may be more efficient data structure
        // arrangements than this (there are!) but
        // this is the algorithm in Wikipedia
        // initialize results array
//        Arrays.fill(sieve, false);


        // the sieve works only for integers > 3, so
        // set these trivially to their proper values
        sieve.set(2);
        sieve.set(3);
//        sieve[0] = false;
//        sieve[1] = false;
//        sieve[2] = true;
//        sieve[3] = true;

        // loop through all possible integer values for x and y
        // up to the square root of the max prime for the sieve
        // we don't need any larger values for x or y since the
        // max value for x or y will be the square root of n
        // in the quadratics
        // the theorem showed that the quadratics will produce all
        // primes that also satisfy their wheel factorizations, so
        // we can produce the value of n from the quadratic first
        // and then filter n through the wheel quadratic
        // there may be more efficient ways to do this, but this
        // is the design in the Wikipedia article
        // loop through all integers for x and y for calculating
        // the quadratics
        for (int x = 1; x <= limitSqrt; x++) {
            for (int y = 1; y <= limitSqrt; y++) {
                // first quadratic using m = 12 and r in R1 = {r : 1, 5}
                int n = (4 * x * x) + (y * y);
                if (n <= limit && (n % 12 == 1 || n % 12 == 5)) {
//                    sieve[n] = !sieve[n];
                    sieve.flip(n);
                }
                // second quadratic using m = 12 and r in R2 = {r : 7}
                n = (3 * x * x) + (y * y);
                if (n <= limit && (n % 12 == 7)) {
//                    sieve[n] = !sieve[n];
                    sieve.flip(n);
                }
                // third quadratic using m = 12 and r in R3 = {r : 11}
                n = (3 * x * x) - (y * y);
                if (x > y && n <= limit && (n % 12 == 11)) {
//                    sieve[n] = !sieve[n];
                    sieve.flip(n);
                } // end if
                // note that R1 union R2 union R3 is the set R
                // R = {r : 1, 5, 7, 11}
                // which is all values 0 < r < 12 where r is
                // a relative prime of 12
                // Thus all primes become candidates
            } // end for
        } // end for
        // remove all perfect squares since the quadratic
        // wheel factorization filter removes only some of them
        for (int n = 5; n <= limitSqrt; n++) {
//            if (sieve[n]) {
            if (sieve.get(n)) {
                int x = n * n;
                for (int i = x; i <= limit; i += x) {
//                    sieve[i] = false;
                    sieve.clear(i);
                } // end for
            } // end if
        } // end for
        // put the results to the System.out device
        // in 10x10 blocks
        for (int i = 0, j = 0; i <= limit; i++) {
//            if (sieve[i]) {
            if (sieve.get(i)) {
                if(!org.apache.commons.math3.primes.Primes.isPrime(i))
                    throw new RuntimeException();
                System.out.printf("%,8d", i);
                if (++j % 10 == 0) {
                    System.out.println();
                } // end if
                if (j % 100 == 0) {
                    System.out.println();
                } // end if
            } // end if
        } // end for
    } // end main
} // end class SieveOfAtkin