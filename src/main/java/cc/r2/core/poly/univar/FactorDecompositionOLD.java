package cc.r2.core.poly.univar;

import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * A holder for polynomial factors.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class FactorDecompositionOLD<T extends lMutablePolynomialAbstract<T>> implements Iterable<T> {
    /** Overall integer factor (polynomial content) */
    long factor;
    /** Factors */
    final List<T> factors;
    /** Either exponents or distinct-degree powers */
    final TIntArrayList exponents;

    public FactorDecompositionOLD() {
        this(new ArrayList<T>(), new TIntArrayList(), 1L);
    }

    FactorDecompositionOLD(List<T> factors, TIntArrayList exponents, long factor) {
        this.factors = factors;
        this.exponents = exponents;
        this.factor = factor;
    }

    @Override
    public Iterator<T> iterator() {
        return factors.iterator();
    }

    T get(int i) { return factors.get(i); }

    int getExponent(int i) { return exponents.get(i); }

    int size() { return factors.size(); }

    boolean isTrivial() {
        return factors.size() == 1;
    }

    protected String toString(boolean infoAsExponents) {
        return toStringFactorization(factors, exponents, factor, infoAsExponents);
    }

    /** Set the content to the specified value */
    FactorDecompositionOLD<T> setNumericFactor(long factor) {
        this.factor = factor;
        return this;
    }

    /** multiply each exponent by a given factor */
    void raiseExponents(long val) {
        for (int i = exponents.size() - 1; i >= 0; --i)
            exponents.set(i, LongArithmetics.safeToInt(exponents.get(i) * val));
    }

    /** canonical form of factor list */
    @SuppressWarnings("unchecked")
    FactorDecompositionOLD<T> canonicalForm() {
        if (factors.size() == 0)
            return this;
        T[] fTmp = factors.toArray((T[]) Array.newInstance(factors.get(0).getClass(), factors.size())); //<- this is ok here, however shitty java generics...
        int[] eTmp = exponents.toArray();
        for (int i = fTmp.length - 1; i >= 0; --i) {
            T poly = fTmp[i];
            if (poly.isMonomial() && eTmp[i] != 1) {
                int degree = poly.degree;
                poly.ensureCapacity(poly.degree * eTmp[i]);
                poly.data[degree * eTmp[i]] = poly.data[degree];
                poly.data[degree] = 0;
                eTmp[i] = 1;
                assert poly.isMonomial();
            }
        }

        ArraysUtil.quickSort(fTmp, eTmp);
        return new FactorDecompositionOLD<T>(new ArrayList<>(Arrays.asList(fTmp)), new TIntArrayList(eTmp), factor);
    }

    /** add another factor */
    void addFactor(T poly, int exponent) {
        if (poly.isConstant()) {
            if (exponent != 1)
                throw new IllegalArgumentException("exponent != 1");
            setNumericFactor(poly.lc());
        } else {
            factors.add(poly);
            exponents.add(exponent);
        }
    }

    /** add another factor */
    void addAll(FactorDecompositionOLD<T> other) {
        factors.addAll(other.factors);
        exponents.addAll(other.exponents);
    }

    @Override
    public String toString() {
        return toString(true);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FactorDecompositionOLD that = (FactorDecompositionOLD) o;
        return factor == that.factor
                && exponents.equals(that.exponents)
                && factors.equals(that.factors);
    }

    @Override
    public int hashCode() {
        int result = (int) (factor^(factor >>> 32));
        for (int i = 0; i < factors.size(); ++i)
            result ^= factors.get(i).hashCode()^exponents.get(i);
        return result;
    }

    /** pretty print for factorization */
    private static String toStringFactorization(List factors, TIntArrayList exponents, long factor, boolean infoAsExponents) {
        if (factors.isEmpty())
            return Long.toString(factor);

        StringBuilder sb = new StringBuilder();
        if (factor != 1) {
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
    T toPolynomial(T factory) {
        return toPolynomial(factory, false);
    }

    /** multiply DDF factors */
    T toPolynomialIgnoringExponents(T factory) {
        return toPolynomial(factory, true);
    }

    private T toPolynomial(T factory, boolean ignoreExps) {
        T r = factory.createOne().multiply(factor);
        for (int i = 0; i < factors.size(); i++) {
            T tmp = ignoreExps ? factors.get(i) : PolynomialArithmetics.polyPow(factors.get(i), exponents.get(i), true);
            r = r.multiply(tmp);
        }
        return r;
    }

    /** decomposition with single numeric factor */
    static <T extends lMutablePolynomialAbstract<T>> FactorDecompositionOLD<T> oneFactor(long factor) {
        return new FactorDecompositionOLD<T>().setNumericFactor(factor);
    }

    /** decomposition with single factor */
    static <T extends lMutablePolynomialAbstract<T>> FactorDecompositionOLD<T> oneFactor(T poly, long factor) {
        FactorDecompositionOLD<T> ts = new FactorDecompositionOLD<>();
        ts.addFactor(poly, 1);
        ts.setNumericFactor(factor);
        return ts;
    }
}
