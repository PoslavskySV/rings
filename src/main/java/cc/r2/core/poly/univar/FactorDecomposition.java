package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A holder for polynomial factors.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class FactorDecomposition<T extends IMutablePolynomial<T>> implements Iterable<T> {
    /** Factors */
    final List<T> factors;
    /** Either exponents or distinct-degree powers */
    final TIntArrayList exponents;

    FactorDecomposition() {
        this(new ArrayList<>(), new TIntArrayList());
    }

    FactorDecomposition(List<T> factors, TIntArrayList exponents) {
        this.factors = factors;
        this.exponents = exponents;
    }

    @Override
    public Iterator<T> iterator() {
        return factors.iterator();
    }

    /** Returns i-th factor */
    public final T get(int i) { return factors.get(i); }

    /** Exponent of i-th factor */
    public final int getExponent(int i) { return exponents.get(i); }

    /** Number of factors */
    public final int size() { return factors.size(); }

    /** Whether there are no non-trivial factors */
    public final boolean isTrivial() { return size() == 1;}

    /** Sum all exponents */
    public final int sumExponents() {
        return exponents.sum();
    }

    /** Multiply each exponent by a given factor */
    public final void raiseExponents(long val) {
        for (int i = exponents.size() - 1; i >= 0; --i)
            exponents.set(i, LongArithmetics.safeToInt(exponents.get(i) * val));
    }

    /** add another factor */
    public abstract void addFactor(T poly, int exponent);

    /** add another factor */
    public final void addAll(FactorDecomposition<T> other) {
        factors.addAll(other.factors);
        exponents.addAll(other.exponents);
    }

    /** Canonical form of factor list */
    @SuppressWarnings("unchecked")
    public final void canonicalForm() {
        if (factors.size() == 0)
            return;
        T[] fTmp = factors.toArray((T[]) Array.newInstance(factors.get(0).getClass(), factors.size())); //<- this is ok here, however shitty java generics...
        int[] eTmp = exponents.toArray();
        for (int i = fTmp.length - 1; i >= 0; --i) {
            T poly = fTmp[i];
            if (poly.isMonomial() && eTmp[i] != 1) {
                poly = PolynomialArithmetics.polyPow(poly, eTmp[i], false);
//                int degree = poly.degree;
//                poly.ensureCapacity(poly.degree * eTmp[i]);
//                poly.data[degree * eTmp[i]] = poly.data[degree];
//                poly.data[degree] = 0;
//                eTmp[i] = 1;
                assert poly.isMonomial();
            }
        }

        ArraysUtil.quickSort(fTmp, eTmp);
        for (int i = 0; i < fTmp.length; i++) {
            factors.set(i, fTmp[i]);
            exponents.set(i, eTmp[i]);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FactorDecomposition<?> that = (FactorDecomposition<?>) o;

        if (factors != null ? !factors.equals(that.factors) : that.factors != null) return false;
        return exponents != null ? exponents.equals(that.exponents) : that.exponents == null;
    }

    @Override
    public int hashCode() {
        int result = factors != null ? factors.hashCode() : 0;
        result = 31 * result + (exponents != null ? exponents.hashCode() : 0);
        return result;
    }

    abstract BigInteger factorAsBigInt();

    @Override
    public String toString() {
        return toString(true);
    }

    protected String toString(boolean infoAsExponents) {
        return toStringFactorization(factors, exponents, factorAsBigInt(), infoAsExponents);
    }

    /** pretty print for factorization */
    private static String toStringFactorization(List factors, TIntArrayList exponents, BigInteger factor, boolean infoAsExponents) {
        if (factors.isEmpty())
            return factor.toString();

        StringBuilder sb = new StringBuilder();
        if (!factor.isOne()) {
            sb.append(factor);
            if (factors.size() > 0)
                sb.append("*");
        }
        for (int i = 0; ; i++) {
            sb.append("(").append(factors.get(i)).append(")");
            if (infoAsExponents && exponents.get(i) != 1)
                sb.append("^").append(exponents.get(i));
            if (i == factors.size() - 1)
                return sb.toString();
            sb.append("*");
        }
    }

    /** multiply factors */
    public final T toPolynomial(T factory) {
        return toPolynomial(factory, false);
    }

    /** multiply DDF factors */
    public final T toPolynomialIgnoringExponents(T factory) {
        return toPolynomial(factory, true);
    }

    abstract T createInitialFactory(T factory);

    private T toPolynomial(T factory, boolean ignoreExps) {
        T r = createInitialFactory(factory);
        for (int i = 0; i < factors.size(); i++) {
            T tmp = ignoreExps ? factors.get(i) : PolynomialArithmetics.polyPow(factors.get(i), exponents.get(i), true);
            r = r.multiply(tmp);
        }
        return r;
    }
}
