package cc.r2.core.poly;

import cc.r2.core.bigint.BigInteger;
import cc.r2.core.bigint.BigIntegerUtil;

/**
 * Abstract domain which holds perfect power decomposition of its cardinality.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class ADomain<E> implements Domain<E> {
    private static final long serialVersionUID = 1L;

    /**
     * if modulus = a^b, a and b are stored in this array
     */
    private final BigInteger[] perfectPowerDecomposition = new BigInteger[2];
    private volatile boolean initialized = false;

    private void checkPerfectPower() {
        // lazy initialization
        if (!initialized) {
            synchronized (perfectPowerDecomposition) {
                if (initialized)
                    return;

                initialized = true;
                if (cardinality() == null) {
                    perfectPowerDecomposition[0] = null;
                    perfectPowerDecomposition[1] = null;
                    return;
                }

                BigInteger[] ipp = BigIntegerUtil.perfectPowerDecomposition(cardinality());
                if (ipp == null) {
                    // not a perfect power
                    perfectPowerDecomposition[0] = cardinality();
                    perfectPowerDecomposition[1] = BigInteger.ONE;
                    return;
                }
                perfectPowerDecomposition[0] = ipp[0];
                perfectPowerDecomposition[1] = ipp[1];

            }
        }
    }

    @Override
    public boolean isPerfectPower() {
        checkPerfectPower();
        return perfectPowerDecomposition[1] != null && !perfectPowerDecomposition[1].isOne();
    }

    @Override
    public BigInteger perfectPowerBase() {
        checkPerfectPower();
        return perfectPowerDecomposition[0];
    }

    @Override
    public BigInteger perfectPowerExponent() {
        checkPerfectPower();
        return perfectPowerDecomposition[1];
    }
}
