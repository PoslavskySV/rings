package cc.r2.core.combinatorics;

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;

import static cc.r2.core.combinatorics.Combinatorics.arrayComparator;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IntDistinctTuplesPortTest {
    @Test
    public void test1() {
        int[] a1 = {1, 1};
        int[] a2 = {1, 1};
        int[] a3 = {1};
        int[][] aa = {a1, a2, a3};
        IntDistinctTuplesPort dcp = new IntDistinctTuplesPort(aa);
        Assert.assertEquals(dcp.take(), null);
    }

    @Test
    public void test2() {
        int[] a1 = {1, 2};
        int[] a2 = {2, 4};
        int[] a3 = {3, 4, 5};
        int[][] aa = {a1, a2, a3};
        IntDistinctTuplesPort dcp = new IntDistinctTuplesPort(aa);
        int[] c;
        Set<int[]> expected = new TreeSet<>(arrayComparator);
        expected.add(new int[]{1, 2, 3});
        expected.add(new int[]{1, 2, 4});
        expected.add(new int[]{1, 2, 5});
        expected.add(new int[]{1, 4, 3});
        expected.add(new int[]{1, 4, 5});
        expected.add(new int[]{2, 4, 3});
        expected.add(new int[]{2, 4, 5});
        Set<int[]> actual = new TreeSet<>(arrayComparator);
        while ((c = dcp.take()) != null)
            actual.add(c.clone());
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void test3() {
        int[] a1 = {1, 2, 3};
        int[] a2 = {1, 2, 3};
        int[] a3 = {1, 2, 3};
        int[][] aa = {a1, a2, a3};
        IntDistinctTuplesPort dcp = new IntDistinctTuplesPort(aa);
        int[] c;
        Set<int[]> actual = new TreeSet<>(arrayComparator);
        while ((c = dcp.take()) != null)
            actual.add(c.clone());

        int[] arr = {1, 2, 3};
        Set<int[]> expected = new TreeSet<>(arrayComparator);
        for (int[] a : Combinatorics.createIntGenerator(3, 3))
            expected.add(permute(arr, a));

        Assert.assertEquals(actual, expected);
    }

    @Test
    public void test4() {
        int[] a1 = {1};
        int[] a2 = {3};
        int[] a3 = {1, 2};
        int[][] aa = {a1, a2, a3};
        IntDistinctTuplesPort dcp = new IntDistinctTuplesPort(aa);
        Assert.assertTrue(Arrays.equals(new int[]{1, 3, 2}, dcp.take()));
        Assert.assertTrue(dcp.take() == null);
    }

    @Test
    public void test5() {
        int[] a1 = {1, 2, 3};
        int[] a2 = {2, 3};
        IntDistinctTuplesPort dcp = new IntDistinctTuplesPort(a1, a2);
        int[] tuple;
        Set<int[]> produced = new TreeSet<>(arrayComparator);
        Set<int[]> expected = new TreeSet<>(arrayComparator);
        expected.add(new int[]{1, 2});
        expected.add(new int[]{1, 3});
        expected.add(new int[]{2, 3});
        expected.add(new int[]{3, 2});
        while ((tuple = dcp.take()) != null)
            produced.add(tuple.clone());
        Assert.assertEquals(expected, produced);
    }

    static int[] permute(int[] array, final int[] permutation) {
        int[] newArray = new int[array.length];
        for (int i = 0; i < permutation.length; ++i)
            newArray[i] = array[permutation[i]];
        return newArray;
    }
}
