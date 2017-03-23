package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import gnu.trove.list.array.TIntArrayList;

import java.util.List;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class bFactorDecomposition<T extends bMutablePolynomialAbstract<T>> extends FactorDecomposition<T> {
    /** overall factor */
    BigInteger factor;

    public bFactorDecomposition(BigInteger factor) {
        super();
        this.factor = factor;
    }

    public bFactorDecomposition(List<T> factors, TIntArrayList exponents, BigInteger factor) {
        super(factors, exponents);
        this.factor = factor;
    }

    /** Set the content to the specified value */
    public bFactorDecomposition<T> setNumericFactor(BigInteger factor) {
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
        return factor;
    }

    /** decomposition with single numeric factor */
    static <T extends bMutablePolynomialAbstract<T>> bFactorDecomposition<T> oneFactor(BigInteger factor) {
        return new bFactorDecomposition<>(factor);
    }

    /** decomposition with single factor */
    static <T extends bMutablePolynomialAbstract<T>> bFactorDecomposition<T> oneFactor(T poly, BigInteger factor) {
        bFactorDecomposition<T> ts = new bFactorDecomposition<>(factor);
        ts.addFactor(poly, 1);
        return ts;
    }
}
