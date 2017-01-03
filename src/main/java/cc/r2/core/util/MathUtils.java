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
package cc.r2.core.util;

import java.util.Arrays;

public final class MathUtils {

    private MathUtils() {
    }

    /**
     * Sort array & return array with removed repetitive values.
     *
     * @param values input array (this method will quickSort this array)
     *
     * @return sorted array of distinct values
     */
    public static int[] getSortedDistinct(int[] values) {
        if (values.length == 0)
            return values;
        Arrays.sort(values);
        int shift = 0;
        int i = 0;
        while (i + shift + 1 < values.length)
            if (values[i + shift] == values[i + shift + 1])
                ++shift;
            else {
                values[i] = values[i + shift];
                ++i;
            }
        values[i] = values[i + shift];
        return Arrays.copyOf(values, i + 1);
    }

    /**
     * Return the set difference B - A for int sets A and B.<br/> Sets A and B
     * must be represented as two sorted int arrays.<br/> Repetitive values in A
     * or B not allowed.
     *
     * @param a sorted array of distinct values. (set A)
     * @param b sorted array of distinct values. (set B)
     *
     * @return the set of elements in B but not in A
     */
    public static int[] intSetDifference(int[] a, int[] b) {
        int bPointer = 0, aPointer = 0;
        int counter = 0;
        while (aPointer < a.length && bPointer < b.length)
            if (a[aPointer] == b[bPointer]) {
                aPointer++;
                bPointer++;
            } else if (a[aPointer] < b[bPointer])
                aPointer++;
            else if (a[aPointer] > b[bPointer]) {
                counter++;
                bPointer++;
            }
        counter += b.length - bPointer;
        int[] result = new int[counter];
        counter = 0;
        aPointer = 0;
        bPointer = 0;
        while (aPointer < a.length && bPointer < b.length)
            if (a[aPointer] == b[bPointer]) {
                aPointer++;
                bPointer++;
            } else if (a[aPointer] < b[bPointer])
                aPointer++;
            else if (a[aPointer] > b[bPointer])
                result[counter++] = b[bPointer++];
        System.arraycopy(b, bPointer, result, counter, b.length - bPointer);
        return result;
    }

    /**
     * Return the union B + A for integer sets A and B.<br/> Sets A and B must
     * be represented as two sorted integer arrays.<br/> Repetitive values in A
     * or B not allowed.
     *
     * @param a sorted array of distinct values. (set A)
     * @param b sorted array of distinct values. (set B)
     *
     * @return the set of elements from B and from A
     */
    public static int[] intSetUnion(int[] a, int[] b) {
        int bPointer = 0, aPointer = 0;
        int counter = 0;
        while (aPointer < a.length && bPointer < b.length)
            if (a[aPointer] == b[bPointer]) {
                aPointer++;
                bPointer++;
                counter++;
            } else if (a[aPointer] < b[bPointer]) {
                aPointer++;
                counter++;
            } else if (a[aPointer] > b[bPointer]) {
                counter++;
                bPointer++;
            }
        counter += (a.length - aPointer) + (b.length - bPointer); //Assert aPoiner==a.length || bPointer==b.length
        int[] result = new int[counter];
        counter = 0;
        aPointer = 0;
        bPointer = 0;
        while (aPointer < a.length && bPointer < b.length)
            if (a[aPointer] == b[bPointer]) {
                result[counter++] = b[bPointer];
                aPointer++;
                bPointer++;
            } else if (a[aPointer] < b[bPointer])
                result[counter++] = a[aPointer++];
            else if (a[aPointer] > b[bPointer])
                result[counter++] = b[bPointer++];
        if (aPointer == a.length)
            System.arraycopy(b, bPointer, result, counter, b.length - bPointer);
        else
            System.arraycopy(a, aPointer, result, counter, a.length - aPointer);
        return result;
    }
}
