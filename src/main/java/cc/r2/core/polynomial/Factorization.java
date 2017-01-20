package cc.r2.core.polynomial;

import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.util.Arrays;
import java.util.List;

/**
 * Polynomial factorization
 */
public final class Factorization {
    /** Integer factor (polynomial content) */
    final long factor;
    /** Factors */
    final MutablePolynomial[] factors;
    /** Either exponents or distinct-degree powers */
    final int[] exponents;

    Factorization(MutablePolynomial[] factors, int[] exponents, long factor) {
        this.factors = factors;
        this.exponents = exponents;
        this.factor = factor;
    }

    Factorization(List<MutablePolynomial> factors, TIntArrayList exponents, long factor) {
        this(factors.toArray(new MutablePolynomial[factors.size()]), exponents.toArray(), factor);
    }

    protected String toString(boolean infoAsExponents) {
        return toStringFactorization(factors, exponents, factor, infoAsExponents);
    }

    Factorization setFactor(long factor) {
        if (factor == this.factor) return this;
        return new Factorization(factors, exponents, factor);
    }

    void raiseExponents(int val) {
        for (int i = exponents.length - 1; i >= 0; --i)
            exponents[i] *= val;
    }

    Factorization canonical() {
        MutablePolynomial[] fTmp = factors.clone();
        int[] eTmp = exponents.clone();
        for (int i = fTmp.length - 1; i >= 0; --i) {
            MutablePolynomial poly = fTmp[i];
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
        return new Factorization(fTmp, eTmp, factor);
    }

    Factorization addFactor(MutablePolynomial poly, int exponent) {
        MutablePolynomial[] factors = Arrays.copyOf(this.factors, this.factors.length + 1);
        int[] exponents = Arrays.copyOf(this.exponents, this.exponents.length + 1);
        factors[factors.length - 1] = poly;
        exponents[exponents.length - 1] = exponent;
        return new Factorization(factors, exponents, factor);
    }

    @Override
    public String toString() {
        return toString(true);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Factorization that = (Factorization) o;

        if (factor != that.factor) return false;
        // Probably incorrect - comparing Object[] arrays with Arrays.equals
        if (!Arrays.equals(factors, that.factors)) return false;
        return Arrays.equals(exponents, that.exponents);
    }

    @Override
    public int hashCode() {
        int result = (int) (factor^(factor >>> 32));
        result = 31 * result + Arrays.hashCode(factors);
        result = 31 * result + Arrays.hashCode(exponents);
        return result;
    }

    /** pretty print for factorization */
    private static String toStringFactorization(Object[] factors, int[] exponents, long factor, boolean infoAsExponents) {
        StringBuilder sb = new StringBuilder();
        if (factor != 1) {
            sb.append(factor);
            if (factors.length > 0)
                sb.append("*");
        }
        for (int i = 0; ; i++) {
            sb.append("(").append(factors[i]).append(")");
            if (infoAsExponents && exponents[i] != 1)
                sb.append("^").append(exponents[i]);
            if (i == factors.length - 1)
                return sb.toString();
            sb.append("*");
        }
    }

    static Factorization oneFactor(long factor) {
        return new Factorization(new MutablePolynomial[0], new int[0], factor);
    }

    static Factorization oneFactor(MutablePolynomial poly, long factor) {
        return new Factorization(new MutablePolynomial[]{poly}, new int[]{1}, factor);
    }
}
