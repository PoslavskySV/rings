package cc.r2.core.combinatorics;

import java.util.Iterator;

/**
 * Parent interface for combinatoric iterators
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class IntCombinatorialGenerator
        implements Iterable<int[]>, Iterator<int[]> {
    /**
     * Resets the iteration
     */
    public abstract void reset();

    /**
     * Returns the reference on the current iteration element.
     *
     * @return the reference on the current iteration element
     */
    public abstract int[] getReference();

    @Override
    public Iterator<int[]> iterator() {
        return this;
    }
}
