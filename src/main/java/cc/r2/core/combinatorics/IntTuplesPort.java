package cc.r2.core.combinatorics;

import java.util.Arrays;

/**
 * This class allows to iterate over all N-tuples (not necessary to be distinct), which can be
 * chosen from {@code N} arrays of integers of the form
 * <i>array</i><sub>i</sub> = [0, 1, 2, ..., K<sub>i</sub>].
 * <br></br>For example, if a set of arrays length is {K<sub>i</sub>} = [2,3,2],
 * then the following tuples will be produced
 * <code><pre>
 *    arr      lastUpdateDepth
 * [0, 0, 0]          0
 * [0, 0, 1]          2
 * [0, 1, 0]          1
 * [0, 1, 1]          2
 * [0, 2, 0]          1
 * [0, 2, 1]          2
 * [1, 0, 0]          0
 * [1, 0, 1]          2
 * [1, 1, 0]          1
 * [1, 1, 1]          2
 * [1, 2, 0]          1
 * [1, 2, 1]          2
 * </pre></code>
 * <p/>
 * <p>This class is implemented via output port pattern and the calculation of the next
 * tuple occurs only on the invocation of {@link #take()}.</p>
 * <p><b>Note:</b> method {@link #take()} returns the same reference on each invocation.
 * So, if it is needed not only to obtain the information from {@link #take()}, but also save the result,
 * it is necessary to clone the returned array.</p>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntTuplesPort implements IntCombinatorialPort {

    private final int[] upperBounds;
    private int[] current;
    private int lastUpdateDepth = -1;

    public IntTuplesPort(final int... upperBounds) {
        checkWithException(upperBounds);
        this.upperBounds = upperBounds;
        this.current = new int[upperBounds.length];
        this.current[upperBounds.length - 1] = -1;
    }

    private static void checkWithException(int[] upperBounds) {
        for (int i : upperBounds)
            if (i < 0)
                throw new IllegalArgumentException("Upper bound cannot be negative.");
    }

    @Override
    public int[] take() {
        int pointer = upperBounds.length - 1;
        boolean next = false;
        ++current[pointer];
        if (current[pointer] == upperBounds[pointer]) {
            current[pointer] = 0;
            next = true;
        }
        while (next && --pointer >= 0) {
            next = false;
            ++current[pointer];
            if (current[pointer] == upperBounds[pointer]) {
                current[pointer] = 0;
                next = true;
            }
        }
        if (lastUpdateDepth != -1)
            lastUpdateDepth = pointer;
        else
            lastUpdateDepth = 0;
        if (next)
            return null;
        return current;
    }

    public int getLastUpdateDepth() {
        return lastUpdateDepth;
    }

    /**
     * Resets the iteration
     */
    public void reset() {
        Arrays.fill(current, 0);
        current[upperBounds.length - 1] = -1;
    }

    @Override
    public int[] getReference() {
        return current;
    }
}
