/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.rings.util;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class ArraysUtilTest {

    @Test
    public void testReverse1() throws Exception {
        long[] arr = {1, 2, 3};
        ArraysUtil.reverse(arr, 0, arr.length);
        assertArrayEquals(new long[]{3, 2, 1}, arr);
        arr = new long[]{1, 2, 3, 4};
        ArraysUtil.reverse(arr, 0, arr.length);
        assertArrayEquals(new long[]{4, 3, 2, 1}, arr);
    }

    @Test
    public void testShort1() {
        short[] target = {2, 1, 0};
        assertArrayEquals(new int[]{2, 1, 0}, ArraysUtil.quickSortP(target));
    }

    @Test
    public void testShort2() {
        short[] target = {2};
        assertArrayEquals(new int[]{0}, ArraysUtil.quickSortP(target));
    }

    @Test
    public void testShort3() {
        short[] target = new short[0];
        assertArrayEquals(new int[0], ArraysUtil.quickSortP(target));
    }

    public static int[] randomPermutation(final int n, RandomGenerator rnd) {
        int[] p = new int[n];
        for (int i = 0; i < n; ++i)
            p[i] = i;
        for (int i = n; i > 1; --i)
            ArraysUtil.swap(p, i - 1, rnd.nextInt(i));
        for (int i = n; i > 1; --i)
            ArraysUtil.swap(p, i - 1, rnd.nextInt(i));
        return p;
    }

    public static int[] permute(int[] array, final int[] permutation) {
        if (array.length != permutation.length)
            throw new IllegalArgumentException();
        int[] newArray = new int[array.length];
        for (int i = 0; i < permutation.length; ++i)
            newArray[i] = array[permutation[i]];
        return newArray;
    }

    @Test
    public void testSortPermutation1() {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100; ++i) {
            int[] a = randomPermutation(10, rnd);
            int[] sorted = a.clone();
            int[] permutation = ArraysUtil.quickSortP(sorted);
            assertArrayEquals(permute(a, permutation), sorted);
        }
    }

    @Test
    public void testBijection1() {
        Integer[] from = {1, 3, 1};
        Integer[] to = {1, 3, 1};
        int[] bijection = {0, 1, 2};
        Assert.assertArrayEquals(bijection, ArraysUtil.bijection(from, to));
    }

    @Test
    public void testBijection2() {
        Integer[] from = {1, 3, 1};
        Integer[] to = {3, 1, 1};
        int[] bijection = {1, 0, 2};
        Assert.assertArrayEquals(bijection, ArraysUtil.bijection(from, to));
    }

    @Test
    public void testQuickSortComparator1() {
        IntComparator comparator = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return Integer.compare(b, a);
            }
        };
        RandomGenerator rnd = new Well1024a();
        int[] array = new int[1000];
        for (int t = 0; t < 100; ++t) {
            for (int i = 0; i < array.length; ++i)
                array[i] = rnd.nextInt(10000);

            ArraysUtil.quickSort(array, comparator);
            for (int i = 1; i < array.length; ++i)
                assertTrue(array[i - 1] >= array[i]);
        }
    }


    @Test
    public void testInsert() throws Exception {
        int[] arr = {0, 1, 2, 3};
        assertArrayEquals(new int[]{99, 0, 1, 2, 3}, ArraysUtil.insert(arr, 0, 99));
        assertArrayEquals(new int[]{0, 99, 1, 2, 3}, ArraysUtil.insert(arr, 1, 99));
        assertArrayEquals(new int[]{0, 1, 2, 3, 99}, ArraysUtil.insert(arr, 4, 99));
    }

//    @Test
//    public void testQuickSortWithCosortAndIntComparator1() {
//        final int degree = 1000;
//        final int[] base = new int[50];
//        for (int i = 0; i < base.length; ++i)
//            base[i] = CC.getRandomGenerator().nextInt(base.length);
//
//        IntComparator comparator = new InducedOrdering(base);
//
//        final int[] array = new int[1000];
//        final int[] cosort = new int[1000];
//
//        for (int i = 1; i < array.length; ++i)
//            cosort[i] = i;
//
//        for (int t = 0; t < 100; ++t) {
//            for (int i = 0; i < array.length; ++i)
//                array[i] = CC.getRandomGenerator().nextInt(degree);
//
//            final int[] _array_ = array.clone(),
//                    _cosort_ = cosort.clone();
//
//            ArraysUtil.quickSort(_array_, _cosort_, comparator);
//
//            for (int i = 0; i < array.length; ++i)
//                assertEquals(array[_cosort_[i]], _array_[i]);
//        }
//    }
}
