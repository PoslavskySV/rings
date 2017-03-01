package cc.r2.core.combinatorics;

/**
 * This class represents an iterator over over all possible unique
 * combinations with permutations (i.e. {0,1} and {1,0} both will appear in the iteration) of {@code k} numbers, which
 * can be chosen from the set of {@code n} numbers (0,1,2,...,{@code n}). The total number of such combinations will be
 * {@code n!/(n-k)!}.
 *
 * <p>For example, for {@code k=2} and {@code n=3}, it will produce
 * the following arrays sequence: [0,1], [1,0], [0,2], [2,0], [1,2],
 * [2,1].</p>
 *
 * <p>The iterator is implemented such that each next combination will be calculated only on
 * the invocation of method {@link #next()} (no pre-calculation of results).</p>
 *
 * <p><b>Note:</b> method {@link #next()} returns the same reference on each invocation.
 * So, if it is needed not only to obtain the information from {@link #next()}, but also to save the result,
 * it is necessary to clone the returned array!</p>
 *
 * <p>Inner implementation of this class is simply uses the combination of {@link IntCombinationsGenerator}
 * and {@link IntPermutationsGenerator}.</p>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @see IntCombinationsGenerator
 * @see IntPermutationsGenerator
 * @since 1.0
 */
public final class IntCombinationPermutationGenerator
        extends IntCombinatorialGenerator
        implements IntCombinatorialPort {
    private final int[] permutation, combination;
    final int[] combinationPermutation;
    private final IntPermutationsGenerator permutationsGenerator;
    private final IntCombinationsGenerator combinationsGenerator;
    private final int k;

    public IntCombinationPermutationGenerator(int n, int k) {
        this.k = k;
        this.combinationsGenerator = new IntCombinationsGenerator(n, k);
        this.combination = this.combinationsGenerator.combination;
        this.permutationsGenerator = new IntPermutationsGenerator(k);
        this.permutation = this.permutationsGenerator.permutation;
        this.combinationPermutation = new int[k];
        combinationsGenerator.next();
        System.arraycopy(combination, 0, combinationPermutation, 0, k);
    }

    @Override
    public int[] take() {
        return hasNext() ? next() : null;
    }

    @Override
    public boolean hasNext() {
        return combinationsGenerator.hasNext() || permutationsGenerator.hasNext();
    }

    /**
     * Calculates and returns the next combination.
     *
     * @return the next combination
     */
    @Override
    public int[] next() {
        if (!permutationsGenerator.hasNext()) {
            permutationsGenerator.reset();
            combinationsGenerator.next();
        }
        permutationsGenerator.next();
        for (int i = 0; i < k; ++i)
            combinationPermutation[i] = combination[permutation[i]];
        return combinationPermutation;
    }

    /**
     * Throws UnsupportedOperationException.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public void reset() {
        permutationsGenerator.reset();
        combinationsGenerator.reset();
        combinationsGenerator.next();
    }

    @Override
    public int[] getReference() {
        return combinationPermutation;
    }
}
