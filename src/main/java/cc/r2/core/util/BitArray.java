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

import java.util.Arrays;

import static java.lang.Math.min;

/**
 * This class represents an "array of booleans" with many fast and useful methods. Consumes ~ 8 times less memory than
 * array of booleans for big sizes. Has slightly different semantics than java's built in {@link java.util.BitSet} and
 * also provides additional functionality like {@link #loadValueFrom(BitArray, int, int, int)} and {@link
 * #copyOfRange(int, int)}.
 */
public final class BitArray {
    public static final BitArray EMPTY = new BitArray(0);

    final int[] data;
    final int size;

    private BitArray(int[] data, int size) {
        this.data = data;
        this.size = size;
    }

    /**
     * Creates an array with specified size. Initial state of all bits is 0 (cleared).
     *
     * @param size
     */
    public BitArray(int size) {
        this.size = size;
        this.data = new int[(size + 31) >>> 5];
    }

    /**
     * Creates a bit array from array of booleans
     *
     * @param array boolean array
     */
    public BitArray(boolean[] array) {
        this(array.length);
        for (int i = 0; i < array.length; ++i)
            if (array[i])
                set(i);
    }

    /**
     * <a href="http://en.wikipedia.org/wiki/Augmented_assignment">Compound assignment</a> and.
     *
     * @param bitArray rhs of and
     */
    public void and(BitArray bitArray) {
        if (bitArray.size != size)
            throw new IllegalArgumentException();

        for (int i = data.length - 1; i >= 0; --i)
            data[i] &= bitArray.data[i];
    }

    /**
     * <a href="http://en.wikipedia.org/wiki/Augmented_assignment">Compound assignment</a> or.
     *
     * @param bitArray rhs of or
     */
    public void or(BitArray bitArray) {
        if (bitArray.size != size)
            throw new IllegalArgumentException();

        for (int i = data.length - 1; i >= 0; --i)
            data[i] |= bitArray.data[i];
    }

    /**
     * <a href="http://en.wikipedia.org/wiki/Augmented_assignment">Compound assignment</a> xor.
     *
     * @param bitArray rhs of xor
     */
    public void xor(BitArray bitArray) {
        if (bitArray.size != size)
            throw new IllegalArgumentException();

        for (int i = data.length - 1; i >= 0; --i)
            data[i] ^= bitArray.data[i];
    }

    /**
     * Inverts all bits in this bit array
     */
    public void not() {
        for (int i = data.length - 1; i >= 0; --i)
            data[i] = ~data[i];

        if (size > 0)
            data[data.length - 1] &= lastElementMask();
    }

    int lastElementMask() {
        if ((size & 0x1F) == 0)
            return 0xFFFFFFFF;
        else
            return 0xFFFFFFFF >>> ((data.length << 5) - size);
        //return ~(0xFFFFFFFF << (32 - ((data.length << 5) - size)));
    }

    /**
     * Returns the number of 1 bits.
     *
     * @return number of 1 bits
     */
    public int bitCount() {
        int bits = 0;
        for (int i : data)
            bits += Integer.bitCount(i);
        return bits;
    }

    /**
     * Clears all bits of this bit array
     */
    public void clearAll() {
        Arrays.fill(data, 0);
    }

    /**
     * Returns a clone of this bit array.
     *
     * @return clone of this bit array
     */
    public BitArray clone() {
        return new BitArray(data.clone(), size);
    }

    //public int[] getBits() {
    //    /*IntArrayList ial = new IntArrayList();
    //    for (int i = 0; i < data.length; ++i)
    //        nextTrailingBit() */
    //    return new int[0];
    //}

    /**
     * Returns {@code true} if there are 1 bits in the same positions. Equivalent to {@code !this.and(other).isEmpty()
     * }
     *
     * @param bitArray other bit array
     * @return {@code true} if there are 1 bits in the same positions
     */
    public boolean intersects(BitArray bitArray) {
        if (bitArray.size != size)
            throw new IllegalArgumentException();

        for (int i = data.length - 1; i >= 0; --i)
            if ((bitArray.data[i] & data[i]) != 0)
                return true;

        return false;
    }

    /**
     * Copy values from the array of the same size
     *
     * @param bitArray bit array to copy values from
     */
    public void loadValueFrom(BitArray bitArray) {
        System.arraycopy(bitArray.data, 0, data, 0, data.length);
    }

    /**
     * Returns the state of specified bit
     *
     * @param i index
     * @return {@code true} if bit is set, {@code false} if bit is cleared
     */
    public boolean get(int i) {
        return (data[i >> 5] & (1 << (i & 0x1F))) != 0;
    }

    /**
     * Sets the specified bit (sets to 1)
     *
     * @param i index of bit
     */
    public void set(int i) {
        data[i >> 5] |= (1 << (i & 0x1F));
    }

    /**
     * Clears the specified bit (sets to 0)
     *
     * @param i index of bit
     */
    public void clear(int i) {
        data[i >> 5] &= ~(1 << (i & 0x1F));
    }

    /**
     * Sets the value of specified bit to specified value
     *
     * @param i     index
     * @param value value
     */
    public void set(int i, boolean value) {
        if (value)
            set(i);
        else
            clear(i);
    }

    /**
     * Sets values at specified positions to specified value
     *
     * @param positions positions
     * @param value     value
     */
    public void setAll(int[] positions, boolean value) {
        for (int i : positions)
            set(i, value);
    }

    /**
     * Sets values at specified positions to specified value
     *
     * @param positions positions
     * @param value     value
     */
    public void setAll(TIntArrayList positions, boolean value) {
        for (int i = positions.size() - 1; i >= 0; --i)
            set(positions.get(i), value);
    }

    /**
     * Sets all bits of this bit array
     */
    public void setAll() {
        for (int i = data.length - 2; i >= 0; --i)
            data[i] = 0xFFFFFFFF;

        if (size > 0)
            data[data.length - 1] = lastElementMask();
    }

    /**
     * Returns the length of this bit array
     *
     * @return length of this bit array
     */
    public int size() {
        return size;
    }

    /**
     * Returns true if all bits in this array are in the 1 state.
     *
     * @return true if all bits in this array are in the 1 state
     */
    public boolean isFull() {
        if (size == 0)
            return true;

        for (int i = data.length - 2; i >= 0; --i)
            if (data[i] != 0xFFFFFFFF)
                return false;

        return data[data.length - 1] == lastElementMask();
    }

    /**
     * Returns true if all bits in this array are in the 0 state.
     *
     * @return true if all bits in this array are in the 0 state
     */
    public boolean isEmpty() {
        for (int i = data.length - 1; i >= 0; --i)
            if (data[i] != 0)
                return false;

        return true;
    }

    /**
     * Returns an array with positions of all "1" bits.
     *
     * @return array with positions of all "1" bits
     */
    public int[] getBits() {
        final int[] bits = new int[bitCount()];
        int i, j = 0;
        if (size < 40 || size >> 2 < bits.length) { // norm ??
            for (i = 0; i < size; ++i)
                if (get(i))
                    bits[j++] = i;
        } else {
            i = -1;
            while ((i = nextBit(i + 1)) != -1) {
                if (!get(i))
                    if (j > 0)
                        nextBit(bits[j - 1] + 1);
                bits[j++] = i;
            }
        }
        assert j == bits.length;
        return bits;
    }

    /**
     * Returns the next "1" bit from the specified position.
     *
     * @param position initial position
     * @return position of the next "1" bit of -1 if all bits after position are 0
     */
    public int nextBit(int position) {
        int ret = position & 0x1F;
        if (ret != 0)
            if ((ret = Integer.numberOfTrailingZeros(data[position >>> 5] >>> ret)) != 32)
                return position + ret;
            else
                position += 32;

        ret = 32;
        position = position >>> 5;
        while (position < data.length &&
                (ret = Integer.numberOfTrailingZeros(data[position++])) == 32) ;

        if (position >= data.length && ret == 32)
            return -1;
        else
            return ((position - 1) << 5) + ret;
    }

    /**
     * Returns the next "0" bit from the specified position (inclusively).
     *
     * @param position initial position
     * @return position of the next "0" bit of -1 if all bits after position are 0
     */
    public int nextZeroBit(int position) {
        //todo review
        int ret = position & 0x1F;
        if (ret != 0)
            if ((ret = Integer.numberOfTrailingZeros((~data[position >>> 5]) >>> ret)) != 32) {
                int r = position + ret;
                return r >= size() ? -1 : r;
            } else
                position += 32;

        ret = 32;
        position = position >>> 5;
        while (position < data.length &&
                (ret = Integer.numberOfTrailingZeros(~data[position++])) == 32) ;

        if (position >= data.length && ret == 32)
            return -1;
        else {
            int r = ((position - 1) << 5) + ret;
            return r >= size() ? -1 : r;
        }
    }

    /**
     * Returns a new bit array, containing values from the specified range
     *
     * @param from lower bound of range
     * @param to   upper bound of range
     * @return new bit array, containing values from the specified range
     */
    public BitArray copyOfRange(int from, int to) {
        BitArray ba = new BitArray(to - from);
        ba.loadValueFrom(this, from, 0, to - from);
        return ba;
    }

    public BitArray copyOfRange(int from) {
        return copyOfRange(from, size);
    }

    /**
     * Analog of {@link System#arraycopy(Object, int, Object, int, int)}, where src is {@code bitArray}.
     *
     * @param bitArray     source
     * @param sourceOffset source offset
     * @param thisOffset   destination offset
     * @param length       number of bits to copy
     */
    public void loadValueFrom(BitArray bitArray, int sourceOffset, int thisOffset, int length) {
        if (bitArray == this)
            throw new IllegalArgumentException("Can't copy from itself.");

        //Forcing alignment on other (bitArray)
        int alignmentOffset = (sourceOffset & 0x1F);

        //Now in words (not in bits)
        sourceOffset = sourceOffset >>> 5;

        if (alignmentOffset != 0) {
            int l = min(32 - alignmentOffset, length);
            loadValueFrom(
                    bitArray.data[sourceOffset] >>> alignmentOffset,
                    thisOffset,
                    l);

            thisOffset += l;
            ++sourceOffset;
            length -= l;
        }

        //Bulk copy
        while (length > 0) {
            loadValueFrom(bitArray.data[sourceOffset], thisOffset, min(32, length));
            length -= 32;
            thisOffset += 32;
            ++sourceOffset;
        }
    }

    /**
     * Load 32 bits or less from single integer
     *
     * @param d        integer with bits
     * @param position offset
     * @param length   length
     */
    void loadValueFrom(int d, int position, int length) {
        if (length == 0)
            return;

        int res = position & 0x1F;
        position = position >>> 5;

        //mask for d
        int mask = 0xFFFFFFFF >>> (32 - length);

        if (res == 0) {
            if (length == 32)
                data[position] = d;
            else {
                data[position] &= ~mask;
                data[position] |= d & mask;
            }
            return;
        }

        data[position] &= ~(mask << res);
        data[position] |= (d & mask) << res;

        length -= (32 - res);
        if (length > 0)
            loadValueFrom(d >>> (32 - res), (position + 1) << 5, length);
    }

    public BitArray append(BitArray other) {
        BitArray ba = new BitArray(size + other.size);
        System.arraycopy(data, 0, ba.data, 0, data.length);
        ba.loadValueFrom(other, 0, size, other.size);
        return ba;
    }

    public BitArray times(int times) {
        BitArray ba = new BitArray(times * size);
        if (times > 0)
            System.arraycopy(data, 0, ba.data, 0, data.length);
        for (int i = 1; i < times; ++i)
            ba.loadValueFrom(this, 0, i * size, size);
        return ba;
    }

    @Override
    public String toString() {
        char[] c = new char[size];
        for (int i = 0; i < size; ++i)
            if (get(i))
                c[i] = '1';
            else
                c[i] = '0';
        return new String(c);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BitArray that = (BitArray) o;

        if (size != that.size) return false;

        return Arrays.equals(data, that.data);
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(data);
        result = 31 * result + size;
        return result;
    }
}
