package cc.r2.core.poly2;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
class GlobalRandom {
    private GlobalRandom() {}

    /** thread local instance of random */
    private static final ThreadLocal<RandomGenerator> ThreadLocalRandom
            = ThreadLocal.withInitial(() -> new Well1024a(0x7f67fcad528cfae9L));

    /** Returns random generator associated with current thread */
    static RandomGenerator getRandom() {
        return ThreadLocalRandom.get();
    }
}
