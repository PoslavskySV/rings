package cc.r2.core.combinatorics;

import gnu.trove.list.array.TIntArrayList;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IntCombinationPermutationGeneratorTest {
    @Test
    public void test1() {
        IntCombinationPermutationGenerator gen = new IntCombinationPermutationGenerator(1, 1);
        Assert.assertTrue(gen.hasNext());
        Assert.assertTrue(Arrays.equals(new int[]{0}, gen.next()));
        Assert.assertTrue(!gen.hasNext());
    }

    @Test
    public void test2() {
        IntCombinationPermutationGenerator gen = new IntCombinationPermutationGenerator(0, 0);
        Assert.assertTrue(gen.hasNext());
        Assert.assertTrue(gen.next().length == 0);
        Assert.assertTrue(!gen.hasNext());
    }

    @Test
    public void test3() {
        TIntArrayList a = new TIntArrayList();
        for (int[] cp : new IntCombinationPermutationGenerator(5, 1)) {
            Assert.assertTrue(cp.length == 1);
            a.add(cp[0]);
        }
        Assert.assertTrue(Arrays.equals(new int[]{0, 1, 2, 3, 4}, a.toArray()));
    }

    @Test
    public void test5() {
        IntCombinationPermutationGenerator gen = new IntCombinationPermutationGenerator(3, 0);
        Assert.assertTrue(gen.hasNext());
        Assert.assertTrue(gen.next().length == 0);
        Assert.assertTrue(!gen.hasNext());
    }
}
