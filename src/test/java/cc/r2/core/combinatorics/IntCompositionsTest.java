package cc.r2.core.combinatorics;

import cc.r2.core.util.ArraysUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.TreeSet;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class IntCompositionsTest {
    @Test
    public void test1() throws Exception {
        for (int integer : new int[]{2, 3, 4, 5, 6})
            for (int nCompositions = 1; nCompositions < 7; nCompositions++) {
                IntCompositions comp = new IntCompositions(integer, nCompositions);
                TreeSet<int[]> set = new TreeSet<>(ArraysUtil.COMPARATOR);
                int count = 0;
                int[] r;
                while ((r = comp.take()) != null) {
                    ++count;
                    set.add(r.clone());
                    Assert.assertEquals(integer, ArraysUtil.sum(r));
                }
                Assert.assertEquals(count, set.size());
            }
    }

    @Test
    public void test2() throws Exception {
        IntCompositions comp = new IntCompositions(3, 4);
        TreeSet<int[]> set = new TreeSet<>(ArraysUtil.COMPARATOR);
        int count = 0;
        int[] r;
        while ((r = comp.take()) != null) {
            ++count;
            set.add(r.clone());
            Assert.assertEquals(3, ArraysUtil.sum(r));
        }
        Assert.assertEquals(count, set.size());
        Assert.assertEquals(20, count);
    }

    @Test
    public void test3() throws Exception {
        IntCompositions comp = new IntCompositions(2, 1);
        TreeSet<int[]> set = new TreeSet<>(ArraysUtil.COMPARATOR);
        int count = 0;
        int[] r;
        while ((r = comp.take()) != null) {
            ++count;
            set.add(r.clone());
            Assert.assertEquals(2, ArraysUtil.sum(r));
        }
        Assert.assertEquals(count, set.size());
        Assert.assertEquals(1, count);
    }

    @Test
    public void test4() throws Exception {
        IntCompositions comp = new IntCompositions(1, 4);
        TreeSet<int[]> set = new TreeSet<>(ArraysUtil.COMPARATOR);
        int count = 0;
        int[] r;
        while ((r = comp.take()) != null) {
            ++count;
            set.add(r.clone());
            Assert.assertEquals(1, ArraysUtil.sum(r));
        }
        Assert.assertEquals(count, set.size());
        Assert.assertEquals(4, count);
    }

    @Test
    public void test5() throws Exception {
        IntCompositions comp = new IntCompositions(2, 4);
        TreeSet<int[]> set = new TreeSet<>(ArraysUtil.COMPARATOR);
        int count = 0;
        int[] r;
        while ((r = comp.take()) != null) {
            ++count;
            set.add(r.clone());
            Assert.assertEquals(2, ArraysUtil.sum(r));
        }
        Assert.assertEquals(count, set.size());
        Assert.assertEquals(10, count);
    }
}