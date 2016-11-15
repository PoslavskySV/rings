package cc.r2.core.number.qsi.factpor;

import gnu.trove.set.hash.TIntHashSet;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by poslavsky on 14/11/2016.
 */
public class PrimesListTest {
    @Test
    public void test1() throws Exception {
        for (int i : PrimesList.SmallPrimes12)
            Assert.assertTrue(org.apache.commons.math3.primes.Primes.isPrime(i));
        for (int i : PrimesList.SmallPrimes10)
            Assert.assertTrue(org.apache.commons.math3.primes.Primes.isPrime(i));

    }

    @Test
    public void test2() throws Exception {
        Assert.assertEquals(PrimesList.SmallPrimes12.length, new TIntHashSet(PrimesList.SmallPrimes12).size());
        Assert.assertEquals(PrimesList.SmallPrimes10.length, new TIntHashSet(PrimesList.SmallPrimes10).size());
    }
}