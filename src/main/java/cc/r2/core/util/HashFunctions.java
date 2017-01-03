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

/**
 * Hashing algorithms. The algorithms was taken from different open sources.
 * <p/>
 * <h3>Links:</h3>
 * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
 * <a href="http://www.burtleburtle.net/bob/hash/doobs.html">http://www.burtleburtle.net/bob/hash/doobs.html</a><br/>
 * <a href="http://bretm.home.comcast.net/~bretm/hash">http://bretm.home.comcast.net/~bretm/hash</a><br/>
 * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
 * <a href="http://en.wikipedia.org/wiki/Jenkins_hash_function's">http://en.wikipedia.org/wiki/Jenkins_hash_function's</a><br/>
 * <a href="http://sites.google.com/site/murmurhash/">http://sites.google.com/site/murmurhash/</a><br/>
 * <a href="http://dmy999.com/article/50/murmurhash-2-java-port">http://dmy999.com/article/50/murmurhash-2-java-port</a><br/>
 * <a href="http://en.wikipedia.org/wiki/MurmurHash">http://en.wikipedia.org/wiki/MurmurHash</a><br/>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class HashFunctions {
    private HashFunctions() {
    }

    /**
     * Robert Jenkins' 96 bit Mix Function.
     * <p/>
     * Variable 'c' contains the input key. When the mixing is complete, variable
     * 'c' also contains the hash result. Variable 'a', and 'b' contain initialized
     * random bits. Notice the total number of internal state is 96 bits, much
     * larger than the final output of 32 bits. Also notice the sequence of
     * subtractions rolls through variable 'a' to variable 'c' three times. Each
     * row will act on one variable, mixing in information from the other two
     * variables, followed by a shift operation.
     * <p/>
     * <p>Subtraction is similar to multiplication in that changes in upper bits
     * of the key do not influence lower bits of the addition. The 9 bit shift
     * operations in Robert Jenkins' mixing algorithm shifts the key to the right
     * 61 bits in total, and shifts the key to the left 34 bits in total. As the
     * calculation is chained, each exclusive-or doubles the number of states.
     * There are at least 2^9 different combined versions of the original key,
     * shifted by various amounts. That is why a single bit change in the key
     * can influence widely apart bits in the hash result.
     * <p/>
     * <p>The uniform distribution of the hash function can be determined from
     * the nature of the subtraction operation. Look at a single bit subtraction
     * operation between a key, and a random bit. If the random bit is 0, then
     * the key remains unchanged. If the random bit is 1, then the key will be
     * flipped. A carry will occur in the case where both the key bit and the
     * random bit are 1. Subtracting the random bits will cause about half of
     * the key bits to be flipped. So even if the key is not uniform, subtracting
     * the random bits will result in uniform distribution.
     * <p/>
     * <h3>Links:</h3>
     * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
     * <a href="http://www.burtleburtle.net/bob/hash/doobs.html">http://www.burtleburtle.net/bob/hash/doobs.html</a><br/>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Jenkins_hash_function's">http://en.wikipedia.org/wiki/Jenkins_hash_function's</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param a initialized random bits
     * @param b initialized random bits
     * @param c key to be hashed
     * @return randomized c bits (hashed c)
     */
    public static int mix(int a, int b, int c) {
        a = a - b;
        a = a - c;
        a = a ^ (c >>> 13);
        b = b - c;
        b = b - a;
        b = b ^ (a << 8);
        c = c - a;
        c = c - b;
        c = c ^ (b >>> 13);
        a = a - b;
        a = a - c;
        a = a ^ (c >>> 12);
        b = b - c;
        b = b - a;
        b = b ^ (a << 16);
        c = c - a;
        c = c - b;
        c = c ^ (b >>> 5);
        a = a - b;
        a = a - c;
        a = a ^ (c >>> 3);
        b = b - c;
        b = b - a;
        b = b ^ (a << 10);
        c = c - a;
        c = c - b;
        c = c ^ (b >>> 15);
        return c;
    }

    /**
     * Based on an original suggestion on Robert Jenkin's part in 1997 and
     * Thomas Wang 2007 updates.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Jenkins_hash_function's">http://en.wikipedia.org/wiki/Jenkins_hash_function's</a><br/>
     *
     * @param key key to be hashed
     * @return hashed value
     */
    public static int JenkinWang32shift(int key) {
        // key = (key << 15) - key - 1;
        // (~x) + y is equivalent to y - x - 1 in two's complement representation.
        key = (~key) + (key << 15);
        key = key ^ (key >>> 12);
        key = key + (key << 2);
        key = key ^ (key >>> 4);
        key = key * 2057; // key = (key + (key << 3)) + (key << 11);
        key = key ^ (key >>> 16);
        return key;
    }

    /**
     * This method uses a combination of bit shifts and integer multiplication
     * to hash the input key.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
     *
     * @param key key to be hashed
     * @return hashed key
     */
    public static int Wang32shiftmult(int key) {
        int c2 = 0x27d4eb2d; // a prime or an odd constant
        key = (key ^ 61) ^ (key >>> 16);
        key = key + (key << 3);
        key = key ^ (key >>> 4);
        key = key * c2;
        key = key ^ (key >>> 15);
        return key;
    }

    /**
     * Based on an original suggestion on Robert Jenkin's part in 1997 and
     * Thomas Wang 2007 updates.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Jenkins_hash_function's">http://en.wikipedia.org/wiki/Jenkins_hash_function's</a><br/>
     *
     * @param key key to be hashed
     * @return hashed value
     */
    public static long JenkinWang64shift(long key) {
        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >>> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >>> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >>> 28);
        key = key + (key << 31);
        return key;
    }

    /**
     * Hashing long to int.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.concentric.net/~ttwang/tech/inthash.htm">http://www.concentric.net/~ttwang/tech/inthash.htm</a><br/>
     *
     * @param key key to be hashed
     * @return hashed value
     */
    public static int Wang64to32shift(long key) {
        key = (~key) + (key << 18); // key = (key << 18) - key - 1;
        key = key ^ (key >>> 31);
        key = key * 21; // key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >>> 11);
        key = key + (key << 6);
        key = key ^ (key >>> 22);
        return (int) key;
    }

    /**
     * Fowler/Noll/Vo hash algorithms FNV_BASIS constant
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param">http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param</a><br/>
     */
    public static final long FNV_BASIS = 0x811c9dc5;
    /**
     * Fowler/Noll/Vo hash algorithms FNV_PRIME constant for 32 bit hash
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param">http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param</a><br/>
     */
    public static final long FNV_PRIME_32 = 16777619;//(1 << 24) + 0x193;
    /**
     * Fowler/Noll/Vo hash algorithms FNV_PRIME constant for 64 bit hash
     * <p/>
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param">http://www.isthe.com/chongo/tech/comp/fnv/#FNV-param</a><br/>
     */
    public static final long FNV_PRIME_64 = 1099511628211L;

    /**
     * Fowler-Noll-Vo 32 bit hash (FNV-1a) for bytes array.<br/>
     * <p/>
     * <h3>Algorithm</h3>
     * <p/>
     * <pre>
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *    hash = hash xor octet_of_data
     *    hash = hash * FNV_prime
     * return hash</pre>
     *
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param bytes bytes array to hash
     * @return hash of the initial bytes array
     */
    public static int FVN32hash(byte[] bytes) {
        long hash = FNV_BASIS;
        for (int i = 0; i < bytes.length; i++) {
            hash ^= 0xFF & bytes[i];
            hash *= FNV_PRIME_32;
        }
        return (int) hash;
    }

    /**
     * Fowler-Noll-Vo 32 bit hash (FNV-1a) for integer key. This is big-endian version (native endianess of JVM).<br/>
     * <p/>
     * <h3>Algorithm</h3>
     * <p/>
     * <pre>
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *    hash = hash xor octet_of_data
     *    hash = hash * FNV_prime
     * return hash</pre>
     *
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param c integer key to be hashed
     * @return hash 32 bit hash
     */
    public static int FVN32hash(int c) {
        long hash = FNV_BASIS;
        hash ^= c >>> 24;
        hash *= FNV_PRIME_32;
        hash ^= 0xFF & (c >>> 16);
        hash *= FNV_PRIME_32;
        hash ^= 0xFF & (c >>> 8);
        hash *= FNV_PRIME_32;
        hash ^= 0xFF & c;
        hash *= FNV_PRIME_32;
        return (int) hash;
    }

    /**
     * Fowler-Noll-Vo 64 bit hash (FNV-1a) for bytes array.<br/>
     * <p/>
     * <h3>Algorithm</h3>
     * <p/>
     * <pre>
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *    hash = hash xor octet_of_data
     *    hash = hash * FNV_prime
     * return hash</pre>
     *
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param bytes bytes array to hash
     * @return hash 64 bit hash
     */
    public static long FVN64hash(byte[] bytes) {
        long hash = FNV_BASIS;
        for (int i = 0; i < bytes.length; i++) {
            hash ^= 0xFF & bytes[i];
            hash *= FNV_PRIME_64;
        }
        return hash;
    }

    /**
     * Fowler-Noll-Vo 64 bit hash (FNV-1a) for long key. This is big-endian version (native endianess of JVM).<br/>
     * <p/>
     * <h3>Algorithm</h3>
     * <p/>
     * <pre>
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *    hash = hash xor octet_of_data
     *    hash = hash * FNV_prime
     * return hash</pre>
     *
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param c long key to be hashed
     * @return hash 64 bit hash
     */
    public static long FVN64hash(long c) {
        long hash = FNV_BASIS;
        hash ^= c >>> 56;
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 48);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 40);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 32);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 24);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 16);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & (c >>> 8);
        hash *= FNV_PRIME_64;
        hash ^= 0xFFL & c;
        hash *= FNV_PRIME_64;
        return hash;
    }

    /**
     * Fowler-Noll-Vo 32 bit hash (FNV-1a) for long key. This is big-endian version (native endianess of JVM).<br/>
     * <p/>
     * <h3>Algorithm</h3>
     * <p/>
     * <pre>
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *    hash = hash xor octet_of_data
     *    hash = hash * FNV_prime
     * return hash</pre>
     *
     * <h3>Links</h3>
     * <a href="http://www.isthe.com/chongo/tech/comp/fnv/">http://www.isthe.com/chongo/tech/comp/fnv/</a><br/>
     * <a href="http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function">http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function</a><br/>
     *
     * @param c long key to be hashed
     * @return hash 32 bit hash
     */
    public static int FVN64to32hash(long c) {
        long hash = FNV_BASIS;
        hash ^= c >>> 56;
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 48);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 40);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 32);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 24);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 16);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & (c >>> 8);
        hash *= FNV_PRIME_32;
        hash ^= 0xFFL & c;
        hash *= FNV_PRIME_32;
        return (int) hash;
    }

    /**
     * MurmurHash hash function for bytes array.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://sites.google.com/site/murmurhash/">http://sites.google.com/site/murmurhash/</a><br/>
     * <a href="http://dmy999.com/article/50/murmurhash-2-java-port">http://dmy999.com/article/50/murmurhash-2-java-port</a><br/>
     * <a href="http://en.wikipedia.org/wiki/MurmurHash">http://en.wikipedia.org/wiki/MurmurHash</a><br/>
     *
     * @param data bytes to be hashed
     * @param seed seed parameter
     * @return 32 bit hash
     */
    @SuppressWarnings("fallthrough")
    public static int MurmurHash2(byte[] data, int seed) {
        // 'm' and 'r' are mixing constants generated offline.
        // They're not really 'magic', they just happen to work well.
        int m = 0x5bd1e995;
        int r = 24;

        // Initialize the hash to a 'random' value
        int len = data.length;
        int h = seed ^ len;

        int i = 0;
        while (len >= 4) {
            int k = data[i] & 0xFF;
            k |= (data[i + 1] & 0xFF) << 8;
            k |= (data[i + 2] & 0xFF) << 16;
            k |= (data[i + 3] & 0xFF) << 24;

            k *= m;
            k ^= k >>> r;
            k *= m;

            h *= m;
            h ^= k;

            i += 4;
            len -= 4;
        }

        switch (len) {
            case 3:
                h ^= (data[i + 2] & 0xFF) << 16;
            case 2:
                h ^= (data[i + 1] & 0xFF) << 8;
            case 1:
                h ^= (data[i] & 0xFF);
                h *= m;
        }

        h ^= h >>> 13;
        h *= m;
        h ^= h >>> 15;

        return h;
    }

    /**
     * MurmurHash hash function integer.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://sites.google.com/site/murmurhash/">http://sites.google.com/site/murmurhash/</a><br/>
     * <a href="http://dmy999.com/article/50/murmurhash-2-java-port">http://dmy999.com/article/50/murmurhash-2-java-port</a><br/>
     * <a href="http://en.wikipedia.org/wiki/MurmurHash">http://en.wikipedia.org/wiki/MurmurHash</a><br/>
     *
     * @param c    int to be hashed
     * @param seed seed parameter
     * @return 32 bit hash
     */
    public static int MurmurHash2(int c, int seed) {
        // 'm' and 'r' are mixing constants generated offline.
        // They're not really 'magic', they just happen to work well.
        int m = 0x5bd1e995;
        // Initialize the hash to a 'random' value
        int h = seed ^ 4;
        c *= m;
        c ^= c >>> 24;
        c *= m;
        h *= m;
        h ^= c;
        h ^= h >>> 13;
        h *= m;
        h ^= h >>> 15;
        return h;
    }

    /**
     * MurmurHash hash function for bytes array with default seed value equals 0x2f1a32b3.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://sites.google.com/site/murmurhash/">http://sites.google.com/site/murmurhash/</a><br/>
     * <a href="http://dmy999.com/article/50/murmurhash-2-java-port">http://dmy999.com/article/50/murmurhash-2-java-port</a><br/>
     * <a href="http://en.wikipedia.org/wiki/MurmurHash">http://en.wikipedia.org/wiki/MurmurHash</a><br/>
     *
     * @param data bytes to be hashed
     * @return 32 bit hash
     */
    public static int MurmurHash2(byte[] data) {
        return MurmurHash2(data, 0x2f1a32b3);
    }

    /**
     * MurmurHash hash function integer with default seed value equals to 0x2f1a32b3.
     * <p/>
     * <h3>Links</h3>
     * <a href="http://sites.google.com/site/murmurhash/">http://sites.google.com/site/murmurhash/</a><br/>
     * <a href="http://dmy999.com/article/50/murmurhash-2-java-port">http://dmy999.com/article/50/murmurhash-2-java-port</a><br/>
     * <a href="http://en.wikipedia.org/wiki/MurmurHash">http://en.wikipedia.org/wiki/MurmurHash</a><br/>
     *
     * @param c int to be hashed\
     * @return 32 bit hash
     */
    public static int MurmurHash2(int c) {
        return MurmurHash2(c, 0x2f1a32b3);
    }
}
