package cc.r2.core.util;

import gnu.trove.set.hash.TIntHashSet;
import junit.framework.Assert;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by poslavsky on 02/01/2017.
 */
public class UtilTest {
    @Test
    public void test1() throws Exception {
        TIntHashSet hashs = new TIntHashSet();
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100000; i++) {
            int[] r = new int[15];
            for (int j = 0; j < r.length; j++) {
                r[j] = rnd.nextInt();
            }
            int hash = ArraysUtil.commutativeHashCode(r);
            Arrays.sort(r);
            Assert.assertEquals(hash, ArraysUtil.commutativeHashCode(r));
            hashs.add(hash);
        }
        System.out.println(hashs.size());
    }
}