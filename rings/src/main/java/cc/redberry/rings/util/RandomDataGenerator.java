package cc.redberry.rings.util;


import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @since 1.0
 */
public class RandomDataGenerator extends org.apache.commons.math3.random.RandomDataGenerator {
    private static final long serialVersionUID = 1L;

    public RandomDataGenerator(RandomGenerator rand) {
        super(rand);
    }

    @Override
    public int nextInt(int lower, int upper) throws NumberIsTooLargeException {
        if (lower == upper)
            return lower;
        return super.nextInt(lower, upper);
    }

    @Override
    public long nextLong(long lower, long upper) throws NumberIsTooLargeException {
        if (lower == upper)
            return lower;
        return super.nextLong(lower, upper);
    }
}
