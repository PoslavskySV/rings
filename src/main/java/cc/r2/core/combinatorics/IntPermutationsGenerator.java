package cc.r2.core.combinatorics;

/**
 * This class represents iterator over all possible permutations of specified
 * dimension written in one-line notation. Number of all permutations of dimension
 * D is D!.
 *
 * <p>Example
 * <code><pre>
 *      IntPermutationsGenerator ig = new IntPermutationsGenerator(3);
 *      while (ig.hasNext())
 *          System.out.println(Arrays.toString(ig.next()))
 * </pre></code>
 * The result will be
 * <code><pre>
 *      [0, 1, 2]
 *      [0, 2, 1]
 *      [1, 0, 2]
 *      [1, 2, 0]
 *      [2, 0, 1]
 *      [2, 1, 0]
 * </pre></code>
 * It is also possible to iterate in the opposite direction via {@link #previous()} method.</p>
 *
 * <p>The iterator is implemented such that each next combination will be calculated only on
 * the invocation of method {@link #next()}.</p>
 *
 * <p><b>Note:</b> method {@link #next()} returns the same reference on each invocation.
 * So, if it is needed not only to obtain the information from {@link #next()}, but also save the result,
 * it is necessary to clone the returned array.</p>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntPermutationsGenerator
        extends IntCombinatorialGenerator
        implements IntCombinatorialPort {

    final int[] permutation;
    private boolean onFirst = true;
    private final int size;

    /**
     * Construct iterator over all permutations with specified dimension starting with identity.
     *
     * @param dimension dimension of permutations
     */
    public IntPermutationsGenerator(int dimension) {
        permutation = new int[dimension];
        for (int i = 0; i < dimension; ++i)
            permutation[i] = i;
        this.size = dimension;
    }

    /**
     * Construct iterator over permutations with specified permutation at
     * the start. If starting permutation is not identity, the iterator will not
     * iterate over all possible permutations, but only from starting permutation up to the
     * last permutation, which is [size-1,size-2,....1,0].
     *
     * <p><b>Note:</b> parameter {@code permutation} is
     * not copied in constructor and the same instance will be used during iteration.</p>
     *
     * @param permutation starting permutation
     * @throws IllegalArgumentException if permutation is inconsistent with
     *                                  <i>one-line</i> notation
     */
    public IntPermutationsGenerator(int[] permutation) {
        this.permutation = permutation;
        this.size = permutation.length;
        for (int i = 0; i < size - 1; ++i) {
            if (permutation[i] >= size || permutation[i] < 0)
                throw new IllegalArgumentException("Wrong permutation input: image of " + i + " element"
                        + " greater then degree");
            for (int j = i + 1; j < size; ++j)
                if (permutation[i] == permutation[j])
                    throw new IllegalArgumentException("Wrong permutation input: to elemets have the same image");
        }
    }

    @Override
    public int[] take() {
        return hasNext() ? next() : null;
    }

    @Override
    public boolean hasNext() {
        return !isLast() || onFirst;
    }

    /**
     * Returns {@code true} if the iteration has more elements, iterating in
     * back order. (In other words, returns {@code true} if {@link #previous()} would
     * return an element rather than throwing an exception.)
     *
     * @return {@code true} if the iteration has more elements
     */
    public boolean hasPrevious() {
        return !isFirst();
    }

    private boolean isLast() {
        for (int i = 0; i < size; i++)
            if (permutation[i] != size - 1 - i)
                return false;
        return true;
    }

    private boolean isFirst() {
        for (int i = 0; i < size; i++)
            if (permutation[i] != i)
                return false;
        return true;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next element in the iteration
     */
    @Override
    public int[] next() {
        if (onFirst) {
            onFirst = false;
            return permutation;
        }
        final int end = size - 1;
        int p = end, low, high, med, s;
        while ((p > 0) && (permutation[p] < permutation[p - 1]))
            p--;
        if (p > 0) //if p==0 then it's the last one
        {
            s = permutation[p - 1];
            if (permutation[end] > s)
                low = end;
            else {
                high = end;
                low = p;
                while (high > low + 1) {
                    med = (high + low) >> 1;
                    if (permutation[med] < s)
                        high = med;
                    else
                        low = med;
                }
            }
            permutation[p - 1] = permutation[low];
            permutation[low] = s;
        }
        high = end;
        while (high > p) {
            med = permutation[high];
            permutation[high] = permutation[p];
            permutation[p] = med;
            p++;
            high--;
        }
        return permutation;
    }

    /**
     * Returns the previous element in the iteration.
     *
     * @return the previous element in the iteration
     */
    public int[] previous() {
        int Nm1 = size - 1;
        int p = Nm1, low, high, s, m;
        while ((p > 0) && (permutation[p] > permutation[p - 1]))
            p--;
        if (p > 0) {
            s = permutation[p - 1];
            if (permutation[Nm1] < s)
                low = Nm1;
            else {
                high = Nm1;
                low = p;
                while (high > low + 1) {
                    m = (high + low) >> 1;
                    if (permutation[m] > s)
                        high = m;
                    else
                        low = m;
                }
            }
            permutation[p - 1] = permutation[low];
            permutation[low] = s;
        }
        high = Nm1;
        while (high > p) {
            m = permutation[high];
            permutation[high] = permutation[p];
            permutation[p] = m;
            p++;
            high--;
        }
        return permutation;
    }

    @Override
    public void reset() {
        onFirst = true;
        for (int i = 0; i < size; ++i)
            permutation[i] = i;
    }

    /**
     * @throws UnsupportedOperationException always
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Returns dimension specified in the constructor
     *
     * @return dimension specified in the constructor
     */
    public int getDimension() {
        return size;
    }

    @Override
    public int[] getReference() {
        return permutation;
    }
}
