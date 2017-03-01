package cc.r2.core.combinatorics;

import java.util.Comparator;

/**
 * This class provides factory and utility methods for combinatorics infrastructure.
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Combinatorics {
    private Combinatorics() {}

    /**
     * <p>Returns an {@link IntCombinatorialGenerator} object, which allows to iterate over
     * all possible unique combinations with permutations (i.e. {0,1} and {1,0} both appears for {@code k=2}) of
     * {@code k} numbers, which can be chosen from the set of {@code n} numbers, numbered in the order
     * 0,1,2,...,{@code n}. The total number of such combinations will be {@code n!/(n-k)!}.</p>
     * <p/>
     * <p>For example, for {@code k=2} and {@code n=3}, this method will produce an iterator over
     * the following arrays: [0,1], [1,0], [0,2], [2,0], [1,2], [2,1].</p>
     *
     * @param n number of elements in the set
     * @param k sample size
     * @return an iterator over all combinations (with permutations) to choose k numbers from n numbers.
     * @see IntCombinatorialGenerator
     */
    public static IntCombinatorialGenerator createIntGenerator(int n, int k) {
        if (n < k)
            throw new IllegalArgumentException();
        if (n == k)
            return new IntPermutationsGenerator(n);
        else
            return new IntCombinationPermutationGenerator(n, k);
    }

    static final Comparator<int[]> arrayComparator = (o1, o2) -> {
        int comp = Integer.compare(o1.length, o2.length);
        if (comp != 0)
            return comp;
        for (int i = 0; i < o1.length; i++)
            if ((comp = Integer.compare(o1[i], o2[i])) != 0)
                return comp;
        return 0;
    };
}
