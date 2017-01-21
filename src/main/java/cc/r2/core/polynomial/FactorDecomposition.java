package cc.r2.core.polynomial;

import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A holder for polynomial factors.
 */
public final class FactorDecomposition {
    /** Integer factor (polynomial content) */
    final long factor;
    /** Factors */
    final List<MutablePolynomial> factors;
    /** Either exponents or distinct-degree powers */
    final TIntArrayList exponents;

    FactorDecomposition(List<MutablePolynomial> factors, TIntArrayList exponents, long factor) {
        this.factors = factors;
        this.exponents = exponents;
        this.factor = factor;
    }

    protected String toString(boolean infoAsExponents) {
        return toStringFactorization(factors, exponents, factor, infoAsExponents);
    }

    FactorDecomposition setFactor(long factor) {
        if (factor == this.factor) return this;
        return new FactorDecomposition(factors, exponents, factor);
    }

    FactorDecomposition raiseExponents(long val) {
        for (int i = exponents.size() - 1; i >= 0; --i)
            exponents.set(i, LongArithmetics.toInt(exponents.get(i) * val));
        return this;
    }

    FactorDecomposition canonical() {
        MutablePolynomial[] fTmp = factors.toArray(new MutablePolynomial[factors.size()]);
        int[] eTmp = exponents.toArray();
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
        return new FactorDecomposition(new ArrayList<>(Arrays.asList(fTmp)), new TIntArrayList(eTmp), factor);
    }

    FactorDecomposition addFactor(MutablePolynomial poly, int exponent) {
        factors.add(poly);
        exponents.add(exponent);
        return this;
    }

    @Override
    public String toString() {
        return toString(true);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FactorDecomposition that = (FactorDecomposition) o;
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

    public MutablePolynomial toPolynomial() {
        return toPolynomial(false);
    }

    MutablePolynomial toPolynomialIgnoringExponents() {
        return toPolynomial(true);
    }

    private MutablePolynomial toPolynomial(boolean ignoreExps) {
        MutablePolynomial r = MutablePolynomial.create(factor);
        for (int i = 0; i < factors.size(); i++) {
            MutablePolynomial tmp = ignoreExps ? factors.get(i) : PolynomialArithmetics.polyPow(factors.get(i), exponents.get(i), true);
            r = r.multiply(tmp);
        }
        return r;
    }

    MutablePolynomial toPolynomial(long modulus) {
        return toPolynomial(modulus, false);
    }

    MutablePolynomial toPolynomialIgnoringExponents(long modulus) {
        return toPolynomial(modulus, true);
    }

    private MutablePolynomial toPolynomial(long modulus, boolean ignoreExps) {
        MutablePolynomial r = MutablePolynomial.create(factor);
        for (int i = 0; i < factors.size(); i++) {
            MutablePolynomial tmp = ignoreExps ? factors.get(i) : PolynomialArithmetics.polyPowMod(factors.get(i), exponents.get(i), modulus, true);
            r = r.multiply(tmp, modulus);
        }
        return r;
    }

    static FactorDecomposition oneFactor(long factor) {
        return new FactorDecomposition(new ArrayList<>(), new TIntArrayList(), factor);
    }

    static FactorDecomposition oneFactor(MutablePolynomial poly, long factor) {
        ArrayList<MutablePolynomial> fs = new ArrayList<>();
        fs.add(poly);
        TIntArrayList exs = new TIntArrayList(new int[]{1});
        return new FactorDecomposition(fs, exs, factor);
    }
}
