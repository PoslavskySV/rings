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

import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.random.BitsStreamGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;

public class BitArrayTest {
    @Test
    public void test1() throws Exception {
        BitArray ba = new BitArray(10);

        assertEquals("0000000000", ba.toString());

        ba.set(3);
        ba.set(5);
        ba.set(8);
        ba.set(5, false);

        assertEquals("0001000010", ba.toString());

        assertEquals(2, ba.bitCount());

        ba.not();
        assertEquals("1110111101", ba.toString());

        assertEquals(8, ba.bitCount());

        ba.and(ba);

        assertEquals("1110111101", ba.toString());

        assertEquals(8, ba.bitCount());

        ba.or(ba);

        assertEquals("1110111101", ba.toString());

        assertEquals(8, ba.bitCount());

        ba.xor(ba);

        assertTrue(ba.isEmpty());

        ba.not();

        assertTrue(ba.isFull());
    }

    @Test
    public void testSetValueFrom1() throws Exception {
        BitArray ba = new BitArray(32);
        ba.set(3);
        ba.set(5);
        ba.set(8);
        ba.set(21);
        ba.set(28);

        int d = ba.data[0];
        String str = ba.toString();

        BitArray arr = new BitArray(128);
        for (int i = 0; i < 64; ++i) {
            arr.loadValueFrom(d, i, 32);
            String expected = chars(i, '0') + str + chars(128 - i - 32, '0');
            assertEquals("On " + i, expected, arr.toString());
        }
    }

    @Test
    public void testSetValueFrom2() throws Exception {
        RandomGenerator rg = new Well19937c(2031);

        for (int k = 0; k < 100; ++k) {
            BitArray ba = new BitArray(32);
            ba.data[0] = rg.nextInt();
            int d = ba.data[0];
            String str = ba.toString();

            BitArray arr = new BitArray(128);
            for (int i = 0; i < 64; ++i) {
                arr.loadValueFrom(d, i, 32);
                String expected = chars(i, ((d & 1) == 1) ? '1' : '0') + str + chars(128 - i - 32, '0');
                assertEquals("On " + i, expected, arr.toString());
            }
        }
    }

    @Test
    public void testSetValueFrom2a() throws Exception {
        RandomGenerator rg = new Well19937c(2031);

        for (int k = 0; k < 100; ++k) {
            BitArray ba = new BitArray(32);
            ba.data[0] = rg.nextInt();
            int d = ba.data[0];
            String str = ba.toString();

            int size = rg.nextInt(512);
            boolean[] arr1 = new boolean[size];
            for (int i = 0; i < size; ++i)
                arr1[i] = rg.nextBoolean();

            BitArray arr = new BitArray(arr1);
            String initial = arr.toString();

            for (int i = 0; i < size - 32; ++i) {
                arr.loadValueFrom(d, i, 32);
                String expected = chars(i, ((d & 1) == 1) ? '1' : '0') + str + initial.substring(32 + i);
                assertEquals("On " + i, expected, arr.toString());
            }
        }
    }

    @Test
    public void testSetValueFrom2b() throws Exception {
        RandomGenerator rg = new Well19937c(2031);

        for (int k = 0; k < 100; ++k) {
            BitArray ba = new BitArray(32);
            ba.data[0] = rg.nextInt();
            int d = ba.data[0];
            String str = ba.toString().substring(0, 16);

            int size = rg.nextInt(512);
            boolean[] arr1 = new boolean[size];
            for (int i = 0; i < size; ++i)
                arr1[i] = rg.nextBoolean();

            BitArray arr = new BitArray(arr1);
            String initial = arr.toString();

            for (int i = 0; i < size - 32; ++i) {
                arr.loadValueFrom(d, i, 16);
                String expected = chars(i, ((d & 1) == 1) ? '1' : '0') + str + initial.substring(16 + i);
                assertEquals("On " + i, expected, arr.toString());
            }
        }
    }

    @Test
    public void testSetValueFrom2c() throws Exception {
        RandomGenerator rg = new Well19937c(2031);

        for (int k = 0; k < 100; ++k) {
            BitArray ba = new BitArray(32);
            ba.data[0] = rg.nextInt();
            int d = ba.data[0];
            String str = ba.toString();

            int size = rg.nextInt(512);
            boolean[] arr1 = new boolean[size];
            for (int i = 0; i < size; ++i)
                arr1[i] = rg.nextBoolean();

            BitArray arr = new BitArray(arr1);
            String initial = arr.toString();

            for (int i = 0; i < size - 32; ++i) {
                arr.loadValueFrom(d, i, 0);
                String expected = initial;
                assertEquals("On " + i, expected, arr.toString());
            }
        }
    }

    @Test
    public void testSetValueFrom3() throws Exception {
        RandomGenerator rg = new Well19937c(2031);

        int offset1, offset2, length;

        for (int k = 0; k < 100; ++k) {
            int size = rg.nextInt(512);

            boolean[] arr1 = new boolean[size], arr2 = new boolean[size];
            for (int i = 0; i < size; ++i)
                arr1[i] = rg.nextBoolean();
            for (int i = 0; i < size; ++i)
                arr2[i] = rg.nextBoolean();

            BitArray ba1 = new BitArray(arr1),
                    ba2 = new BitArray(arr2);


            for (int i = 0; i < 100; ++i) {
                offset1 = rg.nextInt(size);
                offset2 = rg.nextInt(size);
                length = rg.nextInt(size - Math.max(offset1, offset2));
                System.arraycopy(arr1, offset1, arr2, offset2, length);
                ba2.loadValueFrom(ba1, offset1, offset2, length);

                assertTrue(testNormal(ba2));
                assertEquals("On :" + i + ", " + k, new BitArray(arr2), ba2);
            }
        }
    }

    @Test
    public void testCopyOfRange1() throws Exception {
        RandomGenerator rg = new Well19937c(203);

        int offset1, length;

        for (int k = 0; k < 100; ++k) {
            int size = rg.nextInt(512);

            boolean[] arr1 = new boolean[size];
            for (int i = 0; i < size; ++i)
                arr1[i] = rg.nextBoolean();
            BitArray ba1 = new BitArray(arr1), ba2;


            for (int i = 0; i < 100; ++i) {
                offset1 = rg.nextInt(size);
                length = rg.nextInt(size - offset1);
                ba2 = ba1.copyOfRange(offset1, offset1 + length);
                assertTrue(testNormal(ba2));
                assertEquals("On :" + i + ", " + k, new BitArray(Arrays.copyOfRange(arr1, offset1, offset1 + length)), ba2);
            }
        }
    }

    @Test
    public void testNextBit() throws Exception {
        BitArray ba = new BitArray(145);
        ba.set(3);
        ba.set(5);
        ba.set(8);
        ba.set(21);
        ba.set(28);
        ba.set(43);

        assertEquals(3, ba.nextBit(1));
        assertEquals(3, ba.nextBit(3));
        assertEquals(5, ba.nextBit(4));
        assertEquals(28, ba.nextBit(28));
    }

    @Test
    public void testBits1() {
        BitsStreamGenerator random = new Well19937c(325);
        for (int stb = 0; stb < 10000; ++stb) {
            int length;
            boolean[] array = new boolean[length = random.nextInt(200)];
            BitArray bitArray = new BitArray(length);

            int i, bitCount = 0, size;
            TIntArrayList bitsPositions = new TIntArrayList();
            for (i = 0; i < length; ++i)
                if (array[i] = random.nextBoolean()) {
                    bitCount++;
                    bitArray.set(i);
                    bitsPositions.add(i);
                }

            assertEquals(bitCount, bitArray.bitCount());
            assertEquals(bitCount, bitsPositions.size());

            if (bitArray.size() != bitArray.bitCount())
                assertFalse(bitArray.isFull());

            assertArrayEquals(bitsPositions.toArray(), bitArray.getBits());
        }
    }

    @Test
    public void testBits1Sparse() {
        BitsStreamGenerator random = new Well19937c(123);
        for (int stb = 0; stb < 10000; ++stb) {
            int length;
            boolean[] array = new boolean[length = random.nextInt(200)];
            BitArray bitArray = new BitArray(length);

            int i, bitCount = 0, size;
            TIntArrayList bitsPositions = new TIntArrayList();
            for (i = 0; i < length; ++i)
                if (array[i] = (random.nextInt(5) == 0)) {
                    bitCount++;
                    bitArray.set(i);
                    bitsPositions.add(i);
                }

            assertEquals(bitCount, bitArray.bitCount());
            assertEquals(bitCount, bitsPositions.size());

            if (bitArray.size() != bitArray.bitCount())
                assertFalse(bitArray.isFull());

            assertArrayEquals(bitsPositions.toArray(), bitArray.getBits());
        }
    }

    @Test
    public void testPowAppend() throws Exception {
        BitArray ba = new BitArray(10);

        ba.set(3);
        ba.set(5);
        ba.set(8);

        assertEquals("0001010010", ba.toString());

        assertEquals("00010100100001010010000101001000010100100001010010", ba.times(5).toString());

        assertEquals("00010100100001010010", ba.append(ba).toString());
    }

    @Test
    public void testNextZeroBit1() {
        Random r = new Random();
        for (int i = 0; i < 1000; ++i) {
            BitArray bb = randomBitArray(1 + r.nextInt(1000));
            for (int j = 0; j < bb.size(); ++j)
                assertEquals(nextZeroBit(bb, j), bb.nextZeroBit(j));
        }
    }

    @Test
    public void testNextZeroBit1a() {
        BitArray bb = new BitArray(5);
        bb.set(0);
        bb.set(1);
        bb.set(3);
        bb.set(4);
//        System.out.println(bb.nextZeroBit(0));
//        System.out.println(bb.nextZeroBit(1));
//        System.out.println(bb.nextZeroBit(2));
//        System.out.println(bb.nextZeroBit(3));
        System.out.println(bb.nextZeroBit(4));
    }

    private static int nextZeroBit(BitArray array, int position) {
        int i;
        for (i = position; i < array.size(); ++i) {
            if (!array.get(i))
                return i;
        }
        return -1;
    }

    private static BitArray randomBitArray(int size) {
        BitArray b = new BitArray(size);
        Random r = new Random();
        for (int i = (int) (r.nextInt(size) * 0.7); i >= 0; --i)
            b.set(r.nextInt(size));
        return b;
    }

    private boolean testNormal(BitArray ba) {
        if (ba.size == 0)
            return true;
        return (ba.data[ba.data.length - 1] & (~(ba.lastElementMask()))) == 0;
    }

    private String chars(int count, char c) {
        char[] chars = new char[count];
        Arrays.fill(chars, c);
        return new String(chars);
    }
}
