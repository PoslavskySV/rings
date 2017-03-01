package cc.r2.core.combinatorics;

import cc.r2.core.util.OutputPort;

/**
 * This interface is common for all combinatorial iterators.
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface IntCombinatorialPort extends OutputPort<int[]> {
    /**
     * Resets the iteration
     */
    void reset();

    /**
     * Returns the reference to the current iteration element.
     *
     * @return the reference to the current iteration element
     */
    int[] getReference();

    /**
     * Calculates and returns the next combination or null, if no more combinations exist.
     *
     * @return the next combination or null, if no more combinations exist
     */
    @Override
    int[] take();
}
