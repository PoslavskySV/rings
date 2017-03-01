package cc.r2.core.combinatorics;

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IntPermutationsGeneratorTest {

    @Test
    public void test1() {
        class IntArray {
            final int[] data;

            public IntArray(int[] data) {
                this.data = data;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                IntArray intArray = (IntArray) o;

                return Arrays.equals(data, intArray.data);
            }

            @Override
            public int hashCode() {
                return Arrays.hashCode(data);
            }
        }
        IntPermutationsGenerator ig = new IntPermutationsGenerator(8);
        int num = 0;
        Set<IntArray> set = new HashSet<>(40320);
        IntArray a;
        while (ig.hasNext()) {
            ++num;
            a = new IntArray(ig.next().clone());
            Assert.assertTrue(!set.contains(a));
            set.add(a);
        }
        Assert.assertTrue(num == 40320);
    }

    @Test
    public void test2() {
        IntPermutationsGenerator ig = new IntPermutationsGenerator(0);
        Assert.assertTrue(ig.hasNext());
        Assert.assertTrue(ig.next().length == 0);
        Assert.assertTrue(!ig.hasNext());
    }
}
