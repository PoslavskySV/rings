package cc.r2.core.combinatorics;

import java.util.*;

/**
 * This class represents iterator over all possible permutations of specified
 * dimension written in one-line notation (see {@link IntPermutationsGenerator})
 * and allows to specify the niceness of a particular permutations,
 * so they will appear earlier in the iteration if iterator was reset via {@link #reset()}.
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @see IntPermutationsGenerator
 * @since 1.0
 */
public final class IntPriorityPermutationsGenerator implements IntCombinatorialPort {
    private final IntPermutationsGenerator generator;
    private final List<PermutationPriorityTuple> tuples = new ArrayList<>();
    private final Set<PermutationPriorityTuple> set = new HashSet<>();
    private int[] last = null;
    private int lastTuplePointer = 0;

    public IntPriorityPermutationsGenerator(int dimension) {
        generator = new IntPermutationsGenerator(dimension);
    }

    public IntPriorityPermutationsGenerator(int[] initialPermutation) {
        generator = new IntPermutationsGenerator(initialPermutation);
    }

    @Override
    public int[] take() {
        if (lastTuplePointer == tuples.size()) {
            if (!generator.hasNext())
                return null;
            int[] next;
            do {
                if (!generator.hasNext())
                    return null;
                next = generator.next();
            } while (set.contains(new PermutationPriorityTuple(next)));
            last = next;
            return next;
        }
        return tuples.get(lastTuplePointer++).permutation;
    }

    /**
     * Increase niceness of the last returned permutation.
     */
    public void nice() {
        if (last == null) {
            int index = lastTuplePointer - 1;
            int nPriority = ++tuples.get(index).priority;
            int position = index;
            while (--position >= 0 && tuples.get(position).priority < nPriority) ;
            ++position;
            swap(position, index);
            return;
        }
        PermutationPriorityTuple tuple = new PermutationPriorityTuple(last.clone());
        set.add(tuple);
        tuples.add(tuple);
        ++lastTuplePointer;
    }

    @Override
    public void reset() {
        generator.reset();
        lastTuplePointer = 0;
        last = null;
    }

    @Override
    public int[] getReference() {
        return tuples.get(lastTuplePointer - 1).permutation;
    }

    private void swap(int i, int j) {
        PermutationPriorityTuple permutationPriorityTuple = tuples.get(i);
        tuples.set(i, tuples.get(j));
        tuples.set(j, permutationPriorityTuple);
    }

    private static class PermutationPriorityTuple {
        final int[] permutation;
        int priority;

        PermutationPriorityTuple(int[] permutation) {
            this.permutation = permutation;
            this.priority = 1;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            final PermutationPriorityTuple other = (PermutationPriorityTuple) obj;
            return Arrays.equals(this.permutation, other.permutation);
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 89 * hash + Arrays.hashCode(this.permutation);
            return hash;
        }

        @Override
        public String toString() {
            return Arrays.toString(permutation) + " : " + priority;
        }
    }
}
