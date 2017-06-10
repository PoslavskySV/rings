package cc.r2.core.combinatorics;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntCompositions implements IntCombinatorialPort {
    private final int integer;
    private final int nCompositions;
    private final IntCombinationsGenerator generator;
    private int[] array;

    public IntCompositions(int integer, int nCompositions) {
        this.integer = integer;
        this.nCompositions = nCompositions;
        this.generator = new IntCombinationsGenerator(integer + nCompositions - 1, nCompositions - 1);
        this.array = new int[nCompositions];
    }

    @Override
    public void reset() {generator.reset();}

    @Override
    public int[] getReference() {
        return array;
    }

    @Override
    public int[] take() {
        int[] gen = generator.take();
        if (gen == null)
            return array = null;

        if (gen.length == 0) {
            array[0] = integer;
            return array;
        }
        array[0] = gen[0];
        array[array.length - 1] = integer + nCompositions - 1 - gen[gen.length - 1] - 1;
        for (int i = 1; i < array.length - 1; ++i)
            array[i] = gen[i] - gen[i - 1] - 1;

        return array;
    }
}
