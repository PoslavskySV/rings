package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.util.List;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lFactorDecomposition<T extends MutablePolynomialAbstract<T>> extends FactorDecomposition<T> {
    /** overall factor */
    long factor = 1;

    public lFactorDecomposition() {super();}

    public lFactorDecomposition(long factor) {
        super();
        this.factor = factor;
    }

    public lFactorDecomposition(List<T> factors, TIntArrayList exponents, long factor) {
        super(factors, exponents);
        this.factor = factor;
    }

    public lFactorDecomposition(List<T> factors, long factor) {
        this(factors, new TIntArrayList(ArraysUtil.arrayOf(1, factors.size())), factor);
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        lFactorDecomposition<?> that = (lFactorDecomposition<?>) o;

        return factor == that.factor;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (int) (factor^(factor >>> 32));
        return result;
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

//    static bFactorDecomposition<bMutablePolynomialMod> convertZp(lFactorDecomposition<MutablePolynomialMod> decomposition) {
//        bFactorDecomposition<bMutablePolynomialMod> r = new bFactorDecomposition<>();
//        decomposition.factors.forEach(t -> r.factors.add(t.toBigPoly()));
//        r.exponents.addAll(decomposition.exponents);
//        r.setNumericFactor(decomposition.factorAsBigInt());
//        return r;
//    }
//
//    static bFactorDecomposition<bMutablePolynomialZ> convertZ(lFactorDecomposition<MutablePolynomialZ> decomposition) {
//        return convert(decomposition);
//    }

    @SuppressWarnings("unchecked")
    public static <bPoly extends bMutablePolynomialAbstract<bPoly>, lPoly extends MutablePolynomialAbstract<lPoly>>
    bFactorDecomposition<bPoly> convert(lFactorDecomposition<lPoly> decomposition) {
        bFactorDecomposition r = new bFactorDecomposition<bPoly>();
        decomposition.factors.forEach(t -> r.factors.add(t.toBigPoly()));
        r.exponents.addAll(decomposition.exponents);
        r.setNumericFactor(decomposition.factorAsBigInt());
        return r;
    }
}
