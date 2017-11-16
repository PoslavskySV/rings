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


import cc.redberry.rings.bigint.BigInteger;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;
import java.util.function.Function;

/**
 * This class contains additional methods for manipulating arrays (such as sorting and searching). For all quick sort
 * methods the base code was taken from jdk6 {@link Arrays} class.
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @see Arrays
 */
public final class ArraysUtil {
    private ArraysUtil() {
    }

    public static final Comparator<Object> HASH_COMPARATOR = (o1, o2) -> Integer.compare(o1.hashCode(), o2.hashCode());

    /**
     * Lexicographic order
     */
    public static final Comparator<int[]> COMPARATOR = (int[] a, int[] b) -> {
        if (a.length != b.length)
            throw new IllegalArgumentException();
        for (int i = 0; i < a.length; ++i) {
            int c = Integer.compare(a[i], b[i]);
            if (c != 0)
                return c;
        }
        return 0;
    };

    /**
     * Lexicographic order
     */
    public static final Comparator<long[]> COMPARATOR_LONG = (long[] a, long[] b) -> {
        if (a.length != b.length)
            throw new IllegalArgumentException();
        for (int i = 0; i < a.length; ++i) {
            int c = Long.compare(a[i], b[i]);
            if (c != 0)
                return c;
        }
        return 0;
    };

    /**
     * Lexicographic order
     */
    public static final Comparator<Comparable[]> COMPARATOR_GENERIC = (Comparable[] a, Comparable[] b) -> {
        if (a.length != b.length)
            throw new IllegalArgumentException();
        for (int i = 0; i < a.length; ++i) {
            int c = a[i].compareTo(b[i]);
            if (c != 0)
                return c;
        }
        return 0;
    };

    public static int[] flatten(int[][] array) {
        int len = 0;
        for (int[] e : array)
            len += e.length;
        int[] result = new int[len];
        int pointer = 0;
        for (int[] e : array) {
            System.arraycopy(e, 0, result, pointer, e.length);
            pointer += e.length;
        }
        return result;
    }

    public static int[] arrayOf(int val, int len) {
        int[] r = new int[len];
        Arrays.fill(r, val);
        return r;
    }

    public static long[] arrayOf(long val, int len) {
        long[] r = new long[len];
        Arrays.fill(r, val);
        return r;
    }

    public static char[] arrayOf(char val, int len) {
        char[] r = new char[len];
        Arrays.fill(r, val);
        return r;
    }

    public static <T> T[] arrayOf(T val, int len) {
        @SuppressWarnings("unchecked")
        T[] r = (T[]) Array.newInstance(val.getClass(), len);
        Arrays.fill(r, val);
        return r;
    }


    public static int[] negate(int[] arr) {
        for (int i = 0; i < arr.length; i++)
            arr[i] = -arr[i];
        return arr;
    }

    public static long[] negate(long[] arr) {
        for (int i = 0; i < arr.length; i++)
            arr[i] = -arr[i];
        return arr;
    }

    public static BigInteger[] negate(BigInteger[] arr) {
        for (int i = 0; i < arr.length; i++)
            arr[i] = arr[i].negate();
        return arr;
    }

    public static String toString(long[] a, int from, int to) {
        if (a == null)
            return "null";
        int iMax = to - 1;
        if (iMax == -1)
            return "[]";

        StringBuilder b = new StringBuilder();
        b.append('[');
        for (int i = from; ; i++) {
            b.append(a[i]);
            if (i == iMax)
                return b.append(']').toString();
            b.append(", ");
        }
    }

    public static <T> String toString(T[] a, int from, int to) {
        if (a == null)
            return "null";
        int iMax = to - 1;
        if (iMax == -1)
            return "[]";

        StringBuilder b = new StringBuilder();
        b.append('[');
        for (int i = from; ; i++) {
            b.append(a[i]);
            if (i == iMax)
                return b.append(']').toString();
            b.append(", ");
        }
    }

    public static <T> String toString(T[] a, int from, int to, Function<T, String> stringer) {
        if (a == null)
            return "null";
        int iMax = to - 1;
        if (iMax == -1)
            return "[]";

        StringBuilder b = new StringBuilder();
        b.append('[');
        for (int i = from; ; i++) {
            b.append(stringer.apply(a[i]));
            if (i == iMax)
                return b.append(']').toString();
            b.append(", ");
        }
    }

    /**
     * Sort array & return array with removed repetitive values.
     *
     * @param values input array (this method will quickSort this array)
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
     * Sort array & return array with removed repetitive values.
     *
     * @param values input array (this method will quickSort this array)
     * @return sorted array of distinct values
     */
    public static BigInteger[] getSortedDistinct(BigInteger[] values) {
        if (values.length == 0)
            return values;
        Arrays.sort(values);
        int shift = 0;
        int i = 0;
        while (i + shift + 1 < values.length)
            if (values[i + shift].equals(values[i + shift + 1]))
                ++shift;
            else {
                values[i] = values[i + shift];
                ++i;
            }
        values[i] = values[i + shift];
        return Arrays.copyOf(values, i + 1);
    }

    /**
     * Sort array & return array with removed repetitive values.
     *
     * @param values input array (this method will quickSort this array)
     * @return sorted array of distinct values
     */
    public static long[] getSortedDistinct(long[] values) {
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
     * Return the set difference B - A for int sets A and B.<br/> Sets A and B must be represented as two sorted int
     * arrays.<br/> Repetitive values in A or B not allowed.
     *
     * @param main   sorted array of distinct values (set B)
     * @param delete sorted array of distinct values (set A)
     * @return the set of elements in B but not in A
     */
    public static int[] intSetDifference(int[] main, int[] delete) {
        int bPointer = 0, aPointer = 0;
        int counter = 0;
        while (aPointer < delete.length && bPointer < main.length)
            if (delete[aPointer] == main[bPointer]) {
                aPointer++;
                bPointer++;
            } else if (delete[aPointer] < main[bPointer])
                aPointer++;
            else if (delete[aPointer] > main[bPointer]) {
                counter++;
                bPointer++;
            }
        counter += main.length - bPointer;
        int[] result = new int[counter];
        counter = 0;
        aPointer = 0;
        bPointer = 0;
        while (aPointer < delete.length && bPointer < main.length)
            if (delete[aPointer] == main[bPointer]) {
                aPointer++;
                bPointer++;
            } else if (delete[aPointer] < main[bPointer])
                aPointer++;
            else if (delete[aPointer] > main[bPointer])
                result[counter++] = main[bPointer++];
        System.arraycopy(main, bPointer, result, counter, main.length - bPointer);
        return result;
    }

    /**
     * Return the union B + A for integer sets A and B.<br/> Sets A and B must be represented as two sorted integer
     * arrays.<br/> Repetitive values in A or B not allowed.
     *
     * @param a sorted array of distinct values. (set A)
     * @param b sorted array of distinct values. (set B)
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

    public static int[] insert(int[] array, int position, int value) {
        int[] newArray = new int[array.length + 1];
        System.arraycopy(array, 0, newArray, 0, position);
        System.arraycopy(array, position, newArray, position + 1, array.length - position);
        newArray[position] = value;
        return newArray;
    }

    @SuppressWarnings("unchecked")
    public static <T> T[] insert(T[] array, int position, T value) {
        T[] newArray = (T[]) Array.newInstance(value.getClass(), array.length + 1);
        System.arraycopy(array, 0, newArray, 0, position);
        System.arraycopy(array, position, newArray, position + 1, array.length - position);
        newArray[position] = value;
        return newArray;
    }

    public static void reverse(long[] array, int from, int to) {
        for (int i = 0; i < (to - from) / 2; ++i)
            swap(array, from + i, to - 1 - i);
    }

    public static <T> void reverse(T[] array, int from, int to) {
        for (int i = 0; i < (to - from) / 2; ++i)
            swap(array, from + i, to - 1 - i);
    }

    public static <T> void reverse(T[] array) {
        reverse(array, 0, array.length);
    }

    public static void reverse(long[] array) {
        reverse(array, 0, array.length);
    }

    public static int[] short2int(final short[] a) {
        int[] r = new int[a.length];
        for (int i = 0; i < a.length; ++i)
            r[i] = a[i];
        return r;
    }

    public static short[] int2short(final int[] a) {
        short[] r = new short[a.length];
        for (int i = 0; i < a.length; ++i)
            r[i] = (short) a[i];
        return r;
    }


    public static int[] byte2int(final byte[] a) {
        int[] r = new int[a.length];
        for (int i = 0; i < a.length; ++i)
            r[i] = a[i];
        return r;
    }

    public static short[] byte2short(final byte[] a) {
        short[] r = new short[a.length];
        for (int i = 0; i < a.length; ++i)
            r[i] = a[i];
        return r;
    }

    public static byte[] int2byte(final int[] a) {
        byte[] r = new byte[a.length];
        for (int i = 0; i < a.length; ++i)
            r[i] = (byte) a[i];
        return r;
    }

    public static int max(int[] array) {
        int a = Integer.MIN_VALUE;
        for (int i : array)
            a = Math.max(a, i);
        return a;
    }

    public static int max(int[] array, int from , int to) {
        int a = Integer.MIN_VALUE;;
        for (int i = from; i < to; i++)
            a = Math.max(a, array[i]);
        return a;
    }

    public static int[] max(int[] a, int[] b) {
        int[] r = new int[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = Math.max(a[i], b[i]);
        return r;
    }

    public static int min(int[] array) {
        int a = Integer.MAX_VALUE;
        for (int i : array)
            a = Math.min(a, i);
        return a;
    }

    public static int min(int[] array, int from , int to) {
        int a = Integer.MAX_VALUE;
        for (int i = from; i < to; i++)
            a = Math.min(a, array[i]);
        return a;
    }

    public static int[] min(int[] a, int[] b) {
        int[] r = new int[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = Math.min(a[i], b[i]);
        return r;
    }

    public static int firstIndexOf(int element, int[] array) {
        for (int i = 0; i < array.length; i++)
            if (array[i] == element)
                return i;
        return -1;
    }

    public static int firstIndexOf(Object element, Object[] array) {
        for (int i = 0; i < array.length; i++)
            if (element.equals(array[i]))
                return i;
        return -1;
    }

    public static int indexOfMax(int[] array) {
        int index = 0, max = array[index];
        for (int i = 1; i < array.length; i++)
            if (array[i] > max) {
                max = array[i];
                index = i;
            }
        return index;
    }


    public static int[] sequence(int size) {
        return sequence(0, size);
    }

    public static int[] sequence(int from, int to) {
        int[] ret = new int[to - from];
        for (int i = ret.length - 1; i >= 0; --i)
            ret[i] = from + i;
        return ret;
    }

    public static int[][] deepClone(int[][] input) {
        int[][] res = new int[input.length][];
        for (int i = res.length - 1; i >= 0; --i)
            res[i] = input[i].clone();
        return res;
    }

    public static Object[][] deepClone(Object[][] input) {
        Object[][] res = new Object[input.length][];
        for (int i = res.length - 1; i >= 0; --i)
            res[i] = input[i].clone();
        return res;
    }

    public static int sum(final int[] array) {
        return sum(array, 0, array.length);
    }

    public static int sum(final int[] array, int from) {
        return sum(array, from, array.length);
    }

    public static int sum(final int[] array, int from, int to) {
        int s = 0;
        for (int i = from; i < to; ++i)
            s += array[i];
        return s;
    }

    public static int or(final long[] array) {
        return or(array, 0, array.length);
    }

    public static int or(final long[] array, int from) {
        return or(array, from, array.length);
    }

    public static int or(final long[] array, int from, int to) {
        int s = 0;
        for (int i = from; i < to; ++i)
            s |= array[i];
        return s;
    }

    /**
     * This method is similar to {@link #bijection(Comparable[], Comparable[])}  }, but uses specified {@code
     * comparator}.
     *
     * @param from       from array
     * @param to         to array
     * @param comparator comparator
     * @return a bijective mapping from {@code from}-array to {@code to}-array and {@code null} if no mapping exist
     */
    public static <T> int[] bijection(T[] from, T[] to, Comparator<? super T> comparator) {
        //TODO refactor with sorting !!!!
        if (from.length != to.length)
            return null;
        int length = from.length;
        int[] bijection = new int[length];
        Arrays.fill(bijection, -1);
        int i, j;
        OUT:
        for (i = 0; i < length; ++i) {
            for (j = 0; j < length; ++j)
                if (bijection[j] == -1 && comparator.compare(from[i], to[j]) == 0) {
                    bijection[j] = i;
                    continue OUT;
                }
            return null;
        }
        return bijection;
    }

    /**
     * Creates a bijective mapping between two arrays and returns the resulting bijection as array. Method returns null,
     * if no mapping found. <p/>
     * <p>Example: <blockquote><pre>
     *      Integer from[] = {1,2,1,4};
     *      Integer to[] = {2,4,1,1};
     *      int[] bijection = bijection(from,to);
     * </pre></blockquote>
     * <p/> <p> The resulting bijection will be {@code [2, 0, 3, 1]}
     *
     * @param from from array
     * @param to   to array
     * @return a bijective mapping from {@code from}-array to {@code to}-array and {@code null} if no mapping exist
     */
    public static <T extends Comparable<? super T>> int[] bijection(T[] from, T[] to) {
        //TODO refactor with sorting !!!!
        if (from.length != to.length)
            return null;
        int length = from.length;
        int[] bijection = new int[length];
        Arrays.fill(bijection, -1);
        int i, j;
        OUT:
        for (i = 0; i < length; ++i) {
            for (j = 0; j < length; ++j)
                if (bijection[j] == -1 && from[i].compareTo(to[j]) == 0) {
                    bijection[j] = i;
                    continue OUT;
                }
            return null;
        }
        return bijection;
    }

    public static int[] addAll(int[] array1, int... array2) {
        int[] r = new int[array1.length + array2.length];
        System.arraycopy(array1, 0, r, 0, array1.length);
        System.arraycopy(array2, 0, r, array1.length, array2.length);
        return r;
    }

    public static long[] addAll(long[] array1, long... array2) {
        long[] r = new long[array1.length + array2.length];
        System.arraycopy(array1, 0, r, 0, array1.length);
        System.arraycopy(array2, 0, r, array1.length, array2.length);
        return r;
    }

    public static int[] addAll(int[]... arrays) {
        if (arrays.length == 0)
            return new int[0];
        int i, length = 0;
        for (i = 0; i < arrays.length; ++i)
            length += arrays[i].length;
        if (length == 0)
            return new int[0];
        int[] r = new int[length];
        int pointer = 0;
        for (i = 0; i < arrays.length; ++i) {
            System.arraycopy(arrays[i], 0, r, pointer, arrays[i].length);
            pointer += arrays[i].length;
        }
        return r;
    }

    public static int[] remove(int[] array, int i) {
        if (i >= array.length)
            throw new ArrayIndexOutOfBoundsException();
        if (array.length == 1) {
            assert i == 0;
            return new int[0];
        } else if (array.length == 2)
            return new int[]{array[1 ^ i]};
        int[] newArray = new int[array.length - 1];
        System.arraycopy(array, 0, newArray, 0, i);
        if (i != array.length - 1)
            System.arraycopy(array, i + 1, newArray, i, array.length - i - 1);
        return newArray;
    }

    public static long[] remove(long[] array, int i) {
        if (i >= array.length)
            throw new ArrayIndexOutOfBoundsException();
        if (array.length == 1) {
            assert i == 0;
            return new long[0];
        } else if (array.length == 2)
            return new long[]{array[1 ^ i]};
        long[] newArray = new long[array.length - 1];
        System.arraycopy(array, 0, newArray, 0, i);
        if (i != array.length - 1)
            System.arraycopy(array, i + 1, newArray, i, array.length - i - 1);
        return newArray;
    }

    public static <T> T[] remove(T[] array, int i) {
        @SuppressWarnings("unchecked")
        T[] r = (T[]) Array.newInstance(array.getClass().getComponentType(), array.length - 1);
        System.arraycopy(array, 0, r, 0, i);
        if (i < array.length - 1)
            System.arraycopy(array, i + 1, r, i, array.length - i - 1);
        return r;
    }

    /**
     * This code is taken from Apache Commons Lang ArrayUtils. <p/> <p>Adds all the elements of the given arrays into a
     * new array.</p> <p>The new array contains all of the element of {@code array1} followed by all of the elements
     * {@code array2}. When an array is returned, it is always a new array.</p> <p/>
     * <pre>
     * ArrayUtils.addAll(null, null)     = null
     * ArrayUtils.addAll(array1, null)   = cloned copy of array1
     * ArrayUtils.addAll(null, array2)   = cloned copy of array2
     * ArrayUtils.addAll([], [])         = []
     * ArrayUtils.addAll([null], [null]) = [null, null]
     * ArrayUtils.addAll(["a", "b", "c"], ["1", "2", "3"]) = ["a", "b", "c", "1", "2", "3"]
     * </pre>
     *
     * @param <T>    the component type of the array
     * @param array1 the first array whose elements are added to the new array, may be {@code null}
     * @param array2 the second array whose elements are added to the new array, may be {@code null}
     * @return The new array, {@code null} if both arrays are {@code null}. The type of the new array is the type of the
     * first array, unless the first array is null, in which case the type is the same as the second array.
     * @throws IllegalArgumentException if the array types are incompatible
     * @since 2.1
     */
    @SafeVarargs
    public static <T> T[] addAll(T[] array1, T... array2) {
        if (array1 == null)
            return array2.clone();
        else if (array2 == null)
            return array1.clone();
        final Class<?> type1 = array1.getClass().getComponentType();
        @SuppressWarnings("unchecked") // OK, because array is of type T
                T[] joinedArray = (T[]) Array.newInstance(type1, array1.length + array2.length);
        System.arraycopy(array1, 0, joinedArray, 0, array1.length);
        try {
            System.arraycopy(array2, 0, joinedArray, array1.length, array2.length);
        } catch (ArrayStoreException ase) {
            // Check if problem was due to incompatible types
            /*
             * We do this here, rather than before the copy because: - it would
             * be a wasted check most of the time - safer, in case check turns
             * out to be too strict
             */
            final Class<?> type2 = array2.getClass().getComponentType();
            if (!type1.isAssignableFrom(type2))
                throw new IllegalArgumentException("Cannot store " + type2.getName() + " in an array of "
                        + type1.getName(), ase);
            throw ase; // No, so rethrow original
        }
        return joinedArray;
    }

    /**
     * Removes elements at specified {@code positions} in specified {@code array}. This method preserve the relative
     * order of elements in specified {@code array}.
     *
     * @param array     array of elements
     * @param positions positions of elements that should be removed
     * @param <T>       generic type
     * @return new array with removed elements at specified positions
     * @throws ArrayIndexOutOfBoundsException if some position larger then array length
     */
    public static <T> T[] remove(T[] array, int[] positions) {
        if (array == null)
            throw new NullPointerException();
        int[] p = getSortedDistinct(positions);
        if (p.length == 0)
            return array;

        int size = p.length, pointer = 0, s = array.length;
        for (; pointer < size; ++pointer)
            if (p[pointer] >= s)
                throw new ArrayIndexOutOfBoundsException();

        final Class<?> type = array.getClass().getComponentType();
        @SuppressWarnings("unchecked") // OK, because array is of type T
                T[] r = (T[]) Array.newInstance(type, array.length - p.length);

        pointer = 0;
        int i = -1;
        for (int j = 0; j < s; ++j) {
            if (pointer < size - 1 && j > p[pointer])
                ++pointer;
            if (j == p[pointer]) continue;
            else r[++i] = array[j];
        }
        return r;
    }

    /**
     * Removes elements at specified {@code positions} in specified {@code array}. This method preserve the relative
     * order of elements in specified {@code array}.
     *
     * @param array     array of elements
     * @param positions positions of elements that should be removed
     * @return new array with removed elements at specified positions
     * @throws ArrayIndexOutOfBoundsException if some position larger then array length
     */
    public static int[] remove(int[] array, int[] positions) {
        if (array == null)
            throw new NullPointerException();
        int[] p = getSortedDistinct(positions.clone());
        if (p.length == 0)
            return array;

        int size = p.length, pointer = 0, s = array.length;
        for (; pointer < size; ++pointer)
            if (p[pointer] >= s)
                throw new ArrayIndexOutOfBoundsException();

        int[] r = new int[array.length - p.length];

        pointer = 0;
        int i = -1;
        for (int j = 0; j < s; ++j) {
            if (pointer < size - 1 && j > p[pointer])
                ++pointer;
            if (j == p[pointer]) continue;
            else r[++i] = array[j];
        }
        return r;
    }

    /**
     * Selects elements from specified {@code array} at specified {@code positions}. The resulting array preserves the
     * relative order of elements in specified {@code array}.
     *
     * @param array     array of elements
     * @param positions of elements that should be picked out
     * @param <T>       generic type
     * @return the array of elements that picked out from specified positions in specified array
     */
    public static <T> T[] select(T[] array, int[] positions) {
        if (array == null)
            throw new NullPointerException();
        int[] p = getSortedDistinct(positions);
        final Class<?> type = array.getClass().getComponentType();
        @SuppressWarnings("unchecked") // OK, because array is of type T
                T[] r = (T[]) Array.newInstance(type, p.length);
        int i = -1;
        for (int j : p)
            r[++i] = array[j];
        return r;
    }

    /**
     * Converts {@code Set<Integer>} to {@code int[]}
     *
     * @param set a {@link Set} of {@link Integer}
     * @return {@code int[]}
     */
    public static int[] toArray(Set<Integer> set) {
        int i = -1;
        int[] a = new int[set.size()];
        for (Integer ii : set)
            a[++i] = ii;
        return a;
    }

    /**
     * This is the same method to {@link Arrays#binarySearch(int[], int) }. The differs is in the returned value. If key
     * not found, this method returns the position of the first element, witch is closest to key (i.e. if
     * Arrays.binarySearch returns {@code -low-1}, this method returns {@code low}).
     *
     * @param a   the array to be searched
     * @param key the value to be searched for
     * @return index of the search key, if it is contained in the array; otherwise, <tt><i>insertion point</i></tt>. The
     * <i>insertion point</i> is defined as the point at which the key would be inserted into the array: the index of
     * the first element greater than the key, or <tt>a.length</tt> if all elements in the array are less than the
     * specified key.
     */
    public static int binarySearch1(int[] a, int key) {
        return binarySearch1(a, 0, a.length, key);
    }

    /**
     * This is the same method to {@link Arrays#binarySearch(int[], int, int, int) }. The differs is in the returned
     * value. If key not found, this method returns the position of the first element, witch is closest to key (i.e. if
     * Arrays.binarySearch returns {@code -low-1}, this method returns {@code low}).
     *
     * @param a         the array to be searched
     * @param key       the value to be searched for
     * @param fromIndex the index of the first element (inclusive) to be searched
     * @param toIndex   the index of the last element (exclusive) to be searched
     * @return index of the search key, if it is contained in the array; otherwise, <tt><i>insertion point</i></tt>. The
     * <i>insertion point</i> is defined as the point at which the key would be inserted into the array: the index of
     * the first element greater than the key, or <tt>toIndex</tt> if all elements in the array are less than the
     * specified key.
     */
    public static int binarySearch1(int[] a, int fromIndex, int toIndex,
                                    int key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            int midVal = a[mid];

            if (midVal < key)
                low = mid + 1;
            else if (midVal > key)
                high = mid - 1;
            else {
                while (mid > 0 && a[mid - 1] == a[mid]) --mid;
                return mid; // key found
            }
        }
        if (low >= a.length) return low;
        while (low > 0 && a[low - 1] == a[low]) --low;
        return low;  // key not found.
    }

    /**
     * Returns commutative hash code of the data
     *
     * @param data array
     * @return commutative hash
     */
    public static <T> int commutativeHashCode(final T[] data) {
        return commutativeHashCode(data, 0, data.length);
    }

    /**
     * Returns commutative hash code of the data
     *
     * @param data array
     * @return commutative hash
     */
    public static <T> int commutativeHashCode(T[] data, int from, int to) {
        int r = 17;
        for (int i = from; i < to; i++) {
            int h = data[i].hashCode();
            r *= h + 29 ^ h;
        }
        return r;
    }

    /**
     * Returns commutative hash code of the data
     *
     * @param data array
     * @return commutative hash
     */
    public static int commutativeHashCode(final int[] data) {
        return commutativeHashCode(data, 0, data.length);
    }

    /**
     * Returns commutative hash code of the data
     *
     * @param data array
     * @return commutative hash
     */
    public static int commutativeHashCode(int[] data, int from, int to) {
        int r = 17;
        for (int i = from; i < to; i++)
            r *= data[i] + 29 ^ data[i];
        return r;
    }


    /**
     * Sorts the specified array of ints into ascending order using insertion sort algorithm and simultaneously permutes
     * the {@code coSort} ints array in the same way as the target array. This sort guarantee O(n^2) performance in the
     * worst case and O(n) in the best case (nearly sorted input). <p/> <p> This sort is the best choice for small
     * arrays with elements number < 100. <p/> <p>This sort is guaranteed to be <i>stable</i>: equal elements will not
     * be reordered as a result of the sort; <i>adaptive</i>: performance adapts to the initial order of elements and
     * <i>in-place</i>: requires constant amount of additional space.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void insertionSort(int[] target, int[] coSort) {
        insertionSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified array of ints into ascending order using insertion sort algorithm and simultaneously permutes
     * the {@code coSort} ints array in the same way as the target array. This sort guarantee O(n^2) performance in the
     * worst case and O(n) in the best case (nearly sorted input). The range to be sorted extends from index
     * <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the range
     * to be sorted is empty.) <p/> <p> This sort is the best choice for small arrays with elements number < 100. <p/>
     * <p>This sort is guaranteed to be <i>stable</i>: equal elements will not be reordered as a result of the sort;
     * <i>adaptive</i>: performance adapts to the initial order of elements and <i>in-place</i>: requires constant
     * amount of additional space.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void insertionSort(int[] target, int fromIndex, int toIndex, int[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);

        int i, key, j, keyC;
        for (i = fromIndex + 1; i < toIndex; i++) {
            key = target[i];
            keyC = coSort[i];
            for (j = i; j > fromIndex && target[j - 1] > key; j--) {
                target[j] = target[j - 1];
                coSort[j] = coSort[j - 1];
            }
            target[j] = key;
            coSort[j] = keyC;
        }
    }

    /**
     * Sorts the specified array of ints into ascending order using insertion sort algorithm and simultaneously permutes
     * the {@code coSort} longs array in the same way as the specified target array. This sort guarantee O(n^2)
     * performance in the worst case and O(n) in the best case (nearly sorted input). <p/> <p> This sort is the best
     * choice for small arrays with elements number < 100. <p/> <p>This sort is guaranteed to be <i>stable</i>: equal
     * elements will not be reordered as a result of the sort; <i>adaptive</i>: performance adapts to the initial order
     * of elements and <i>in-place</i>: requires constant amount of additional space.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void insertionSort(int[] target, long[] coSort) {
        insertionSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified array of ints into ascending order using insertion sort algorithm and simultaneously permutes
     * the {@code coSort} ints array in the same way as the target array. This sort guarantee O(n^2) performance in the
     * worst case and O(n) in the best case (nearly sorted input). The range to be sorted extends from index
     * <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the range
     * to be sorted is empty.) <p/> <p> This sort is the best choice for small arrays with elements number < 100. <p/>
     * <p>This sort is guaranteed to be <i>stable</i>: equal elements will not be reordered as a result of the sort;
     * <i>adaptive</i>: performance adapts to the initial order of elements and <i>in-place</i>: requires constant
     * amount of additional space.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the specified target array, during sorting
     *                  procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void insertionSort(int[] target, int fromIndex, int toIndex, long[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);

        int i, key, j;
        long keyC;
        for (i = fromIndex + 1; i < toIndex; i++) {
            key = target[i];
            keyC = coSort[i];
            for (j = i; j > fromIndex && target[j - 1] > key; j--) {
                target[j] = target[j - 1];
                coSort[j] = coSort[j - 1];
            }
            target[j] = key;
            coSort[j] = keyC;
        }
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements using insertion sort algorithm and simultaneously permutes the {@code coSort} objects array in the same
     * way then specified target array. This sort guarantee O(n^2) performance in the worst case and O(n) in the best
     * case (nearly sorted input). <p/> <p> This sort is the best choice for small arrays with elements number < 100.
     * <p/> <p>This sort is guaranteed to be <i>stable</i>: equal elements will not be reordered as a result of the
     * sort; <i>adaptive</i>: performance adapts to the initial order of elements and <i>in-place</i>: requires constant
     * amount of additional space.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static <T extends Comparable<T>> void insertionSort(T[] target, Object[] coSort) {
        insertionSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements using insertion sort algorithm and simultaneously permutes the {@code coSort} objects array in the same
     * way then specified target array. This sort guarantee O(n^2) performance in the worst case and O(n) in the best
     * case (nearly sorted input). The range to be sorted extends from index <tt>fromIndex</tt>, inclusive, to index
     * <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the range to be sorted is empty.) <p/> <p> This
     * sort is the best choice for small arrays with elements number < 100. <p/> <p>This sort is guaranteed to be
     * <i>stable</i>: equal elements will not be reordered as a result of the sort; <i>adaptive</i>: performance adapts
     * to the initial order of elements and <i>in-place</i>: requires constant amount of additional space.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static <T extends Comparable<T>> void insertionSort(T[] target, int fromIndex, int toIndex, Object[] coSort) {
        int i, j;
        T key;
        Object keyC;
        for (i = fromIndex + 1; i < toIndex; i++) {
            key = target[i];
            keyC = coSort[i];
            for (j = i; j > fromIndex && target[j - 1].compareTo(key) > 0; j--) {
                target[j] = target[j - 1];
                coSort[j] = coSort[j - 1];
            }
            target[j] = key;
            coSort[j] = keyC;
        }
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements using insertion sort algorithm and simultaneously permutes the {@code coSort} objects array in the same
     * way then specified target array. This sort guarantee O(n^2) performance in the worst case and O(n) in the best
     * case (nearly sorted input). <p/> <p> This sort is the best choice for small arrays with elements number < 100.
     * <p/> <p>This sort is guaranteed to be <i>stable</i>: equal elements will not be reordered as a result of the
     * sort; <i>adaptive</i>: performance adapts to the initial order of elements and <i>in-place</i>: requires constant
     * amount of additional space.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static <T extends Comparable<T>> void insertionSort(T[] target, int[] coSort) {
        insertionSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements using insertion sort algorithm and simultaneously permutes the {@code coSort} objects array in the same
     * way then specified target array. This sort guarantee O(n^2) performance in the worst case and O(n) in the best
     * case (nearly sorted input). The range to be sorted extends from index <tt>fromIndex</tt>, inclusive, to index
     * <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the range to be sorted is empty.) <p/> <p> This
     * sort is the best choice for small arrays with elements number < 100. <p/> <p>This sort is guaranteed to be
     * <i>stable</i>: equal elements will not be reordered as a result of the sort; <i>adaptive</i>: performance adapts
     * to the initial order of elements and <i>in-place</i>: requires constant amount of additional space.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static <T extends Comparable<T>> void insertionSort(T[] target, int fromIndex, int toIndex, int[] coSort) {
        int i, j;
        T key;
        int keyC;
        for (i = fromIndex + 1; i < toIndex; i++) {
            key = target[i];
            keyC = coSort[i];
            for (j = i; j > fromIndex && target[j - 1].compareTo(key) > 0; j--) {
                target[j] = target[j - 1];
                coSort[j] = coSort[j - 1];
            }
            target[j] = key;
            coSort[j] = keyC;
        }
    }

    /**
     * Sorts the specified array of ints into ascending order using TimSort algorithm and simultaneously permutes the
     * {@code coSort} ints array in the same way as the target array. <p/> <p> NOTE: using of this method is very good
     * for large arrays with more then 100 elements, in other case using of insertion sort is highly recommended. <p/>
     * <p>This sort is guaranteed to be <i>stable</i>: equal elements will not be reordered as a result of the sort.
     * <p/> <p> The code was taken from {@link Arrays#sort(java.lang.Object[]) } and adapted for integers. For more
     * information look there.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws ClassCastException if the array contains elements that are not <i>mutually comparable</i> (for example,
     *                            strings and integers)
     * @see Arrays#sort(java.lang.Object[])
     */
    public static void timSort(int target[], int[] coSort) {
        IntTimSort.sort(target, coSort);
    }

    /**
     * Sorts the specified array of ints into ascending order using stable sort algorithm and simultaneously permutes
     * the {@code coSort} ints array in the same way as the target array. If length of specified array is less than 100
     * - insertion sort algorithm performed, otherwise - TimSort. <p/> <p>This sort is guaranteed to be <i>stable</i>:
     * equal elements will not be reordered as a result of the sort.
     *
     * @param target the array to be sorted
     * @param cosort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws ClassCastException       if the array contains elements that are not <i>mutually comparable</i> (for
     *                                  example, strings and integers)
     * @throws IllegalArgumentException if coSort length less then target length.
     * @throws IllegalArgumentException if target == coSort (as references).
     * @see #insertionSort(int[], int[])
     * @see #timSort(int[], int[])
     */
    public static void stableSort(int target[], int[] cosort) {
        if (target.length > 100)
            ArraysUtil.timSort(target, cosort);
        else
            ArraysUtil.insertionSort(target, cosort);
    }

    /**
     * Sorts the specified array and returns the resulting permutation
     *
     * @param target int array
     * @return sorting permutation
     */
    public static int[] quickSortP(int[] target) {
        int[] permutation = new int[target.length];
        for (int i = 1; i < target.length; ++i)
            permutation[i] = i;
        quickSort(target, 0, target.length, permutation);
        return permutation;
    }

    // =================  QUICKSORT INT[] INT[] =================

    /**
     * Sorts the specified target array of ints into ascending numerical order and simultaneously permutes the {@code
     * coSort} ints array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays class. <p/>
     * The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a
     * Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers
     * n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/>
     * <p><b>NOTE: remember this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can
     * be perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like
     * an insertion sort or Tim sort.</b> <p/> <p><b>NOTE:</b> The method throws {@code IllegalArgumentException} if
     * {@code target == coSort}, because in this case no sorting will be perfomed.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static void quickSort(int[] target, int[] coSort) {
        quickSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified range of the specified target array of ints into ascending numerical order and simultaneously
     * permutes the {@code coSort} ints array in the same way as the target array. The range to be sorted extends from
     * index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the
     * range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm
     * is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function",
     * Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n)
     * performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p><b>NOTE:
     * remember this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b> <p/> <p><b>NOTE:</b> The method throws {@code IllegalArgumentException} if {@code
     * target == coSort}, because in this case no sorting will be perfomed.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     * @throws IllegalArgumentException       if target == coSort (as references).
     */
    public static void quickSort(int[] target, int fromIndex, int toIndex, int[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, int[]) }, but without range checking and toIndex ->
     * length (see params). Throws {@code IllegalArgumentException} if {@code target == coSort}, because in this case no
     * sorting will be perfomed . <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the
     * {@code coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use
     * stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static void quickSort1(int target[], int fromIndex, int length, int[] coSort) {
        if (target == coSort)
            throw new IllegalArgumentException("Target reference == coSort reference.");
        quickSort2(target, fromIndex, length, coSort);
    }

    private static void quickSort2(int target[], int fromIndex, int length, int[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1] > target[j]; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b] <= v) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c] >= v) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort2(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort2(target, n - s, s, coSort);

    }

    private static void swap(int x[], int a, int b, int[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    /**
     * Swaps x[a] with x[b].
     */
    public static void swap(int x[], int a, int b) {
        int t = x[a];
        x[a] = x[b];
        x[b] = t;
    }

    private static void vecswap(int x[], int a, int b, int n, int[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    /**
     * Returns the index of the median of the three indexed integers.
     */
    private static int med3(int x[], int a, int b, int c) {
        return (x[a] < x[b]
                ? (x[b] < x[c] ? b : x[a] < x[c] ? c : a)
                : (x[b] > x[c] ? b : x[a] > x[c] ? c : a));
    }

    // =================  QUICKSORT LONG[] LONG[] =================

    /**
     * Sorts the specified target array of ints into ascending numerical order and simultaneously permutes the {@code
     * coSort} longs array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays class. <p/>
     * The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a
     * Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers
     * n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(long[] target, long[] coSort) {
        quickSort1(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified range of the specified target array of ints into ascending numerical order and simultaneously
     * permutes the {@code coSort} longs array in the same way as the target array. The range to be sorted extends from
     * index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the
     * range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm
     * is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function",
     * Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n)
     * performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * performed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(long[] target, int fromIndex, int toIndex, long[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, long[])  ) }, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * performed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     */
    public static void quickSort1(long target[], int fromIndex, int length, long[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1] > target[j]; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        long v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b] <= v) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c] >= v) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, coSort);

    }

    private static void swap(long x[], int a, int b, long[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    private static void vecswap(long x[], int a, int b, int n, long[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    private static int med3(long x[], int a, int b, int c) {
        return (x[a] < x[b]
                ? (x[b] < x[c] ? b : x[a] < x[c] ? c : a)
                : (x[b] > x[c] ? b : x[a] > x[c] ? c : a));
    }

    // =================  QUICKSORT INT[] LONG[] =================

    /**
     * Sorts the specified target array of ints into ascending numerical order and simultaneously permutes the {@code
     * coSort} longs array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays class. <p/>
     * The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a
     * Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers
     * n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(int[] target, long[] coSort) {
        quickSort1(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified range of the specified target array of ints into ascending numerical order and simultaneously
     * permutes the {@code coSort} longs array in the same way as the target array. The range to be sorted extends from
     * index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the
     * range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm
     * is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function",
     * Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n)
     * performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * performed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(int[] target, int fromIndex, int toIndex, long[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, long[])  ) }, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * performed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     */
    public static void quickSort1(int target[], int fromIndex, int length, long[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1] > target[j]; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b] <= v) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c] >= v) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, coSort);

    }

    private static void swap(int x[], int a, int b, long[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    /**
     * Swaps x[a] with x[b].
     */
    public static void swap(long x[], int a, int b) {
        long t = x[a];
        x[a] = x[b];
        x[b] = t;
    }

    private static void vecswap(int x[], int a, int b, int n, long[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    // =================  QUICKSORT OBJECT[] OBJECT[] =================

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements and simultaneously permutes the {@code coSort} objects array in the same way then specified target
     * array. <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm is a tuned quicksort,
     * adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function", Software-Practice and
     * Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n) performance on many data
     * sets that cause other quicksorts to degrade to quadratic performance. <p/> <p><b>NOTE: this is unstable sort
     * algorithm, so additional combinatorics of the {@code coSort} array can be perfomed. Use this method only if you
     * are sure, in what you are doing. If not - use stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if <tt>fromIndex &gt; toIndex</tt>
     * @throws IllegalArgumentException if coSort length less then target length.
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort(T[] target, Object[] coSort) {
        quickSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements and simultaneously permutes the {@code coSort} objects array in the same way then specified target
     * array. The range to be sorted extends from index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>,
     * exclusive. (If <tt>fromIndex==toIndex</tt>, the range to be sorted is empty.)<p> <p/> The code was taken from the
     * jdk6 Arrays class. <p/> The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas
     * McIlroy's "Engineering a Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November
     * 1993). This algorithm offers n*log(n) performance on many data sets that cause other quicksorts to degrade to
     * quadratic performance. <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the
     * {@code coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use
     * stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     * @throws IllegalArgumentException       if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort(T[] target, int fromIndex, int toIndex, Object[] coSort) {
        if (target == coSort)
            throw new IllegalArgumentException();
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(Comparable[], int, int, Object[])}, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort1(T[] target, int fromIndex, int length, Object[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1].compareTo(target[j]) > 0; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        T v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b].compareTo(v) <= 0) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c].compareTo(v) >= 0) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, coSort);

    }

    private static void swap(Object[] x, int a, int b, Object[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    private static void vecswap(Object[] x, int a, int b, int n, Object[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    // =================  QUICKSORT OBJECT[] INT[] =================

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements and simultaneously permutes the {@code coSort} objects array in the same way then specified target
     * array. <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm is a tuned quicksort,
     * adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function", Software-Practice and
     * Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n) performance on many data
     * sets that cause other quicksorts to degrade to quadratic performance. <p/> <p><b>NOTE: this is unstable sort
     * algorithm, so additional combinatorics of the {@code coSort} array can be perfomed. Use this method only if you
     * are sure, in what you are doing. If not - use stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if <tt>fromIndex &gt; toIndex</tt>
     * @throws IllegalArgumentException if coSort length less then target length.
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort(T[] target, int[] coSort) {
        quickSort(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified target array of objects into ascending order, according to the natural ordering of its
     * elements and simultaneously permutes the {@code coSort} objects array in the same way then specified target
     * array. The range to be sorted extends from index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>,
     * exclusive. (If <tt>fromIndex==toIndex</tt>, the range to be sorted is empty.)<p> <p/> The code was taken from the
     * jdk6 Arrays class. <p/> The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas
     * McIlroy's "Engineering a Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November
     * 1993). This algorithm offers n*log(n) performance on many data sets that cause other quicksorts to degrade to
     * quadratic performance. <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the
     * {@code coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use
     * stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     * @throws IllegalArgumentException       if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort(T[] target, int fromIndex, int toIndex, int[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(Comparable[], int, int, Object[])}, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException if target == coSort (as references).
     */
    public static <T extends Comparable<T>> void quickSort1(T[] target, int fromIndex, int length, int[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1].compareTo(target[j]) > 0; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        T v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b].compareTo(v) <= 0) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c].compareTo(v) >= 0) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, coSort);

    }

    private static void swap(Object[] x, int a, int b, int[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    private static void vecswap(Object[] x, int a, int b, int n, int[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    /**
     * Swaps x[a] with x[b].
     */
    public static void swap(Object[] x, int a, int b) {
        Object t = x[a];
        x[a] = x[b];
        x[b] = t;
    }

    private static <T extends Comparable<T>> int med3(T[] x, int a, int b, int c) {
        return (x[a].compareTo(x[b]) < 0
                ? (x[b].compareTo(x[c]) < 0 ? b : x[a].compareTo(x[c]) < 0 ? c : a)
                : (x[b].compareTo(x[c]) > 0 ? b : x[a].compareTo(x[c]) > 0 ? c : a));
    }

    // =================  QUICKSORT INT[] OBJECT[] =================

    /**
     * Sorts the specified target array of ints into ascending numerical order and simultaneously permutes the {@code
     * coSort} Objects array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays class.
     * <p/> The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's
     * "Engineering a Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This
     * algorithm offers n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic
     * performance. <p/> <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code
     * coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use stable
     * sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(int[] target, Object[] coSort) {
        quickSort1(target, 0, target.length, coSort);
    }

    /**
     * Sorts the specified range of the specified target array of ints into ascending numerical order and simultaneously
     * permutes the {@code coSort} Objects array in the same way as the target array. The range to be sorted extends
     * from index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>,
     * the range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting
     * algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort
     * Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers
     * n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(int[] target, int fromIndex, int toIndex, Object[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, Object[])  ) }, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     */
    public static void quickSort1(int target[], int fromIndex, int length, Object[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1] > target[j]; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b] <= v) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c] >= v) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, coSort);

    }

    private static void swap(int x[], int a, int b, Object[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    private static void vecswap(int x[], int a, int b, int n, Object[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    public static int[] quickSortP(short[] target) {
        int[] permutation = new int[target.length];
        for (int i = 1; i < target.length; ++i)
            permutation[i] = i;
        quickSort(target, 0, target.length, permutation);
        return permutation;
    }

    // =================  QUICKSORT SHORT[] INT[] =================

    /**
     * Sorts the specified range of the specified target array of ints into ascending numerical order and simultaneously
     * permutes the {@code coSort} ints array in the same way as the target array. The range to be sorted extends from
     * index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If <tt>fromIndex==toIndex</tt>, the
     * range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays class. <p/> The sorting algorithm
     * is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a Sort Function",
     * Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers n*log(n)
     * performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/> <p><b>NOTE:
     * remember this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b> <p/> <p><b>NOTE:</b> The method throws {@code IllegalArgumentException} if {@code
     * target == coSort}, because in this case no sorting will be perfomed.
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param toIndex   the index of the last element (exclusive) to be sorted
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(short[] target, int fromIndex, int toIndex, int[] coSort) {
        rangeCheck(target.length, fromIndex, toIndex);
        rangeCheck(coSort.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, coSort);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, int[]) }, but without range checking and toIndex ->
     * length (see params). Throws {@code IllegalArgumentException} if {@code target == coSort}, because in this case no
     * sorting will be perfomed . <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the
     * {@code coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use
     * stable sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target    the array to be sorted
     * @param fromIndex the index of the first element (inclusive) to be sorted
     * @param length    the length of the sorting subarray.
     * @param coSort    the array which will be permuted in the same way as the target array, during sorting procedure
     */
    public static void quickSort1(short target[], int fromIndex, int length, int[] coSort) {
        quickSort2(target, fromIndex, length, coSort);
    }

    private static void quickSort2(short target[], int fromIndex, int length, int[] coSort) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && target[j - 1] > target[j]; j--)
                    swap(target, j, j - 1, coSort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s);
                m = med3(target, m - s, m, m + s);
                n = med3(target, n - 2 * s, n - s, n);
            }
            m = med3(target, l, m, n); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && target[b] <= v) {
                if (target[b] == v)
                    swap(target, a++, b, coSort);
                b++;
            }
            while (c >= b && target[c] >= v) {
                if (target[c] == v)
                    swap(target, c, d--, coSort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, coSort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, coSort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, coSort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort2(target, fromIndex, s, coSort);
        if ((s = d - c) > 1)
            quickSort2(target, n - s, s, coSort);

    }

    /**
     * Sorts the specified target array of shorts into ascending numerical order and simultaneously permutes the {@code
     * coSort} ints array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays class. <p/>
     * The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's "Engineering a
     * Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This algorithm offers
     * n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic performance. <p/>
     * <p><b>NOTE: remember this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can
     * be perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like
     * an insertion sort or Tim sort.</b> <p/> <p><b>NOTE:</b> The method throws {@code IllegalArgumentException} if
     * {@code target == coSort}, because in this case no sorting will be perfomed.
     *
     * @param target the array to be sorted
     * @param coSort the array which will be permuted in the same way as the target array during sorting procedure
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(short[] target, int[] coSort) {
        quickSort(target, 0, target.length, coSort);
    }

    private static void swap(short x[], int a, int b, int[] coSort) {
        swap(x, a, b);
        swap(coSort, a, b);
    }

    /**
     * Swaps x[a] with x[b].
     */
    private static void swap(short x[], int a, int b) {
        short t = x[a];
        x[a] = x[b];
        x[b] = t;
    }

    private static void vecswap(short x[], int a, int b, int n, int[] coSort) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b, coSort);
    }

    /**
     * Returns the index of the median of the three indexed integers.
     */
    private static int med3(short x[], int a, int b, int c) {
        return (x[a] < x[b]
                ? (x[b] < x[c] ? b : x[a] < x[c] ? c : a)
                : (x[b] > x[c] ? b : x[a] > x[c] ? c : a));
    }

    ////////////////////////////////////// COMPARATOR /////////////////////////////////////////////////////

    /**
     * Sorts the specified range of the specified target array of ints into order specified by {@link IntComparator}
     * using quicksort.
     *
     * @param target     the array to be sorted
     * @param comparator custom comparator
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(int[] target, IntComparator comparator) {
        quickSort1(target, 0, target.length, comparator);
    }

    /**
     * Sorts the specified range of the specified target array of ints into order specified by {@link IntComparator}
     * using quicksort.
     *
     * @param target     the array to be sorted
     * @param fromIndex  the index of the first element (inclusive) to be sorted
     * @param toIndex    the index of the last element (exclusive) to be sorted
     * @param comparator comparator
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(int[] target, int fromIndex, int toIndex, IntComparator comparator) {
        rangeCheck(target.length, fromIndex, toIndex);
        quickSort1(target, fromIndex, toIndex - fromIndex, comparator);
    }

    /**
     * Sorts the specified range of the specified target array of ints into order specified by {@link IntComparator}
     * using quicksort.
     *
     * @param target     the array to be sorted
     * @param fromIndex  the index of the first element (inclusive) to be sorted
     * @param length     the length of the sorting subarray.
     * @param comparator comparator
     */
    public static void quickSort1(int target[], int fromIndex, int length, IntComparator comparator) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && comparator.compare(target[j - 1], target[j]) > 0; j--)
                    swap(target, j, j - 1);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s, comparator);
                m = med3(target, m - s, m, m + s, comparator);
                n = med3(target, n - 2 * s, n - s, n, comparator);
            }
            m = med3(target, l, m, n, comparator); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && comparator.compare(target[b], v) <= 0) {
                if (comparator.compare(target[b], v) == 0)
                    swap(target, a++, b);
                b++;
            }
            while (c >= b && comparator.compare(target[c], v) >= 0) {
                if (comparator.compare(target[c], v) == 0)
                    swap(target, c, d--);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, comparator);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, comparator);
    }


    /////////////////////////////// QUICK SORT INTCOMPARATOR COSORT ////////////////////////////////////////

    /**
     * Sorts the specified target array of ints according to {@link IntComparator} and simultaneously permutes the
     * {@code coSort} Objects array in the same way as the target array. <p/> The code was taken from the jdk6 Arrays
     * class. <p/> The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's
     * "Engineering a Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This
     * algorithm offers n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic
     * performance. <p/> <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code
     * coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use stable
     * sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target     the array to be sorted
     * @param comparator custom comparator
     * @throws IllegalArgumentException if coSort length less then target length.
     */
    public static void quickSort(int[] target, int[] cosort, IntComparator comparator) {
        quickSort1(target, 0, target.length, cosort, comparator);
    }

    /**
     * Sorts the specified range of the specified target array of ints according to {@link IntComparator} and
     * simultaneously permutes the {@code coSort} Objects array in the same way as the target array. The range to be
     * sorted extends from index <tt>fromIndex</tt>, inclusive, to index <tt>toIndex</tt>, exclusive. (If
     * <tt>fromIndex==toIndex</tt>, the range to be sorted is empty.)<p> <p/> The code was taken from the jdk6 Arrays
     * class. <p/> The sorting algorithm is a tuned quicksort, adapted from Jon L. Bentley and M. Douglas McIlroy's
     * "Engineering a Sort Function", Software-Practice and Experience, Vol. 23(11) P. 1249-1265 (November 1993). This
     * algorithm offers n*log(n) performance on many data sets that cause other quicksorts to degrade to quadratic
     * performance. <p/> <p/> <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code
     * coSort} array can be perfomed. Use this method only if you are sure, in what you are doing. If not - use stable
     * sort methods like an insertion sort or Tim sort.</b>
     *
     * @param target     the array to be sorted
     * @param fromIndex  the index of the first element (inclusive) to be sorted
     * @param toIndex    the index of the last element (exclusive) to be sorted
     * @param comparator comparator
     * @throws IllegalArgumentException       if <tt>fromIndex &gt; toIndex</tt>
     * @throws ArrayIndexOutOfBoundsException if <tt>fromIndex &lt; 0</tt> or <tt>toIndex &gt; target.length</tt> or
     *                                        <tt>toIndex &gt; coSort.length</tt>
     */
    public static void quickSort(int[] target, int fromIndex, int toIndex, int[] cosort, IntComparator comparator) {
        rangeCheck(target.length, fromIndex, toIndex);
        if (target == cosort)
            throw new IllegalArgumentException("Same array references.");
        quickSort1(target, fromIndex, toIndex - fromIndex, cosort, comparator);
    }

    /**
     * This method is the same as {@link #quickSort(int[], int, int, Object[])  ) }, but without range checking. <p/>
     * <p><b>NOTE: this is unstable sort algorithm, so additional combinatorics of the {@code coSort} array can be
     * perfomed. Use this method only if you are sure, in what you are doing. If not - use stable sort methods like an
     * insertion sort or Tim sort.</b>
     *
     * @param target     the array to be sorted
     * @param fromIndex  the index of the first element (inclusive) to be sorted
     * @param length     the length of the sorting subarray.
     * @param comparator comparator
     */
    private static void quickSort1(int target[], int fromIndex, int length, int[] cosort, IntComparator comparator) {
        // Insertion quickSort on smallest arrays
        if (length < 7) {
            for (int i = fromIndex; i < length + fromIndex; i++)
                for (int j = i; j > fromIndex && comparator.compare(target[j - 1], target[j]) > 0; j--)
                    swap(target, j, j - 1, cosort);
            return;
        }

        // Choose a partition element, v
        int m = fromIndex + (length >> 1);       // Small arrays, middle element
        if (length > 7) {
            int l = fromIndex;
            int n = fromIndex + length - 1;
            if (length > 40) {        // Big arrays, pseudomedian of 9
                int s = length / 8;
                l = med3(target, l, l + s, l + 2 * s, comparator);
                m = med3(target, m - s, m, m + s, comparator);
                n = med3(target, n - 2 * s, n - s, n, comparator);
            }
            m = med3(target, l, m, n, comparator); // Mid-size, med of 3
        }
        int v = target[m];

        // Establish Invariant: v* (<v)* (>v)* v*
        int a = fromIndex, b = a, c = fromIndex + length - 1, d = c;
        while (true) {
            while (b <= c && comparator.compare(target[b], v) <= 0) {
                if (comparator.compare(target[b], v) == 0)
                    swap(target, a++, b, cosort);
                b++;
            }
            while (c >= b && comparator.compare(target[c], v) >= 0) {
                if (comparator.compare(target[c], v) == 0)
                    swap(target, c, d--, cosort);
                c--;
            }
            if (b > c)
                break;
            swap(target, b++, c--, cosort);
        }

        // Swap partition elements back to middle
        int s, n = fromIndex + length;
        s = Math.min(a - fromIndex, b - a);
        vecswap(target, fromIndex, b - s, s, cosort);
        s = Math.min(d - c, n - d - 1);
        vecswap(target, b, n - s, s, cosort);

        // Recursively quickSort non-partition-elements
        if ((s = b - a) > 1)
            quickSort1(target, fromIndex, s, cosort, comparator);
        if ((s = d - c) > 1)
            quickSort1(target, n - s, s, cosort, comparator);
    }

    /**
     * Returns the index of the median of the three indexed integers.
     */
    private static int med3(int x[], int a, int b, int c, IntComparator comparator) {
        return (comparator.compare(x[a], x[b]) < 0
                ? (comparator.compare(x[b], x[c]) < 0 ? b : comparator.compare(x[a], x[c]) < 0 ? c : a)
                : (comparator.compare(x[b], x[c]) > 0 ? b : comparator.compare(x[a], x[c]) > 0 ? c : a));
    }

    private static void vecswap(int x[], int a, int b, int n) {
        for (int i = 0; i < n; i++, a++, b++)
            swap(x, a, b);
    }


    ////////////////////////////////////// UTILS ///////////////////////////////////////////////////////////

    /**
     * Check that fromIndex and toIndex are in range, and throw an appropriate exception if they aren't.
     */
    private static void rangeCheck(int arrayLen, int fromIndex, int toIndex) {
        if (fromIndex > toIndex)
            throw new IllegalArgumentException("fromIndex(" + fromIndex
                    + ") > toIndex(" + toIndex + ")");
        if (fromIndex < 0)
            throw new ArrayIndexOutOfBoundsException(fromIndex);
        if (toIndex > arrayLen)
            throw new ArrayIndexOutOfBoundsException(toIndex);
    }
}
