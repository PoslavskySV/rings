package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import gnu.trove.list.array.TIntArrayList;

import java.util.List;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lFactorDecomposition<T extends MutablePolynomialAbstract<T>> extends FactorDecomposition<T> {
    /** overall factor */
    long factor;

    public lFactorDecomposition() {super();}

    public lFactorDecomposition(long factor) {
        super();
        this.factor = factor;
    }

    public lFactorDecomposition(List<T> factors, TIntArrayList exponents, long factor) {
        super(factors, exponents);
        this.factor = factor;
    }

    /** Set the content to the specified value */
    public lFactorDecomposition<T> setNumericFactor(long factor) {
        this.factor = factor;
        return this;
    }

    @Override
    public void addFactor(T poly, int exponent) {
        if (poly.isConstant()) {
            if (exponent != 1)
                throw new IllegalArgumentException("exponent != 1");
            setNumericFactor(poly.lc());
        } else {
            factors.add(poly);
            exponents.add(exponent);
        }
    }

    @Override
    T createInitialFactory(T factory) {
        return factory.createConstant(factor);
    }

    @Override
    BigInteger factorAsBigInt() {
        return BigInteger.valueOf(factor);
    }


    /** decomposition with single numeric factor */
    static <T extends MutablePolynomialAbstract<T>> lFactorDecomposition<T> oneFactor(long factor) {
        return new lFactorDecomposition<>(factor);
    }

    /** decomposition with single factor */
    static <T extends MutablePolynomialAbstract<T>> lFactorDecomposition<T> oneFactor(T poly, long factor) {
        lFactorDecomposition<T> ts = new lFactorDecomposition<>(factor);
        ts.addFactor(poly, 1);
        return ts;
    }
}
