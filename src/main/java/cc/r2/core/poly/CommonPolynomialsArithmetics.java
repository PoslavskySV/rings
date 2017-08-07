package cc.r2.core.poly;

import cc.r2.core.poly.multivar.AMultivariatePolynomial;
import cc.r2.core.poly.multivar.MultivariateGCD;
import cc.r2.core.poly.multivar.MultivariateReduction;
import cc.r2.core.poly.univar.DivisionWithRemainder;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariateGCD;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.stream.StreamSupport;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class CommonPolynomialsArithmetics {
    private CommonPolynomialsArithmetics() {}

    @SuppressWarnings("unchecked")
    public static <Poly extends IGeneralPolynomial<Poly>>
    Poly PolynomialGCD(Poly a, Poly b) {
        if (a instanceof IUnivariatePolynomial)
            return (Poly) UnivariateGCD.PolynomialGCD((IUnivariatePolynomial) a, (IUnivariatePolynomial) b);
        else if (a instanceof AMultivariatePolynomial)
            return (Poly) MultivariateGCD.PolynomialGCD((AMultivariatePolynomial) a, (AMultivariatePolynomial) b);
        else
            throw new RuntimeException();
    }

    @SuppressWarnings("unchecked")
    public static <Poly extends IGeneralPolynomial<Poly>>
    Poly[] PolynomialDivideAndRemainder(Poly a, Poly b) {
        if (a instanceof IUnivariatePolynomial)
            return (Poly[]) DivisionWithRemainder.divideAndRemainder((IUnivariatePolynomial) a, (IUnivariatePolynomial) b, true);
        else if (a instanceof AMultivariatePolynomial)
            return (Poly[]) MultivariateReduction.divideAndRemainder((AMultivariatePolynomial) a, (AMultivariatePolynomial) b);
        else
            throw new RuntimeException();
    }

    public static <Poly extends IGeneralPolynomial<Poly>>
    boolean coprimeQ(Poly... polynomials) {
        for (int i = 0; i < polynomials.length - 1; i++)
            for (int j = i + 1; j < polynomials.length; j++)
                if (!PolynomialGCD(polynomials[i], polynomials[j]).isConstant())
                    return false;
        return true;
    }

    public static <Poly extends IGeneralPolynomial<Poly>>
    boolean coprimeQ(Iterable<Poly> polynomials) {
        if (!polynomials.iterator().hasNext())
            throw new IllegalArgumentException();
        Poly factory = polynomials.iterator().next();
        return coprimeQ(StreamSupport.stream(polynomials.spliterator(), false).toArray(factory::arrayNewInstance));
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e}
     */
    public static <T extends IGeneralPolynomial<T>> T polyPow(final T base, long exponent, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 1 || base.isOne())
            return copy ? base.clone() : base;
        T result = base.createOne();
        T k2p = copy ? base.clone() : base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = result.multiply(k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p);
        }
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @param cache    cache to store all intermediate powers
     * @return {@code base} in a power of {@code e}
     */
    public static <T extends IGeneralPolynomial<T>> T polyPow(final T base, int exponent, boolean copy,
                                                              TIntObjectHashMap<T> cache) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 1)
            return copy ? base.clone() : base;

        T cached = cache.get(exponent);
        if (cached != null)
            return cached.clone();

        T result = base.createOne();
        T k2p = copy ? base.clone() : base;
        int rExp = 0, kExp = 1;
        for (; ; ) {
            if ((exponent&1) != 0)
                cache.put(rExp += kExp, result.multiply(k2p).clone());
            exponent = exponent >> 1;
            if (exponent == 0) {
                cache.put(rExp, result);
                return result;
            }
            cache.put(kExp *= 2, k2p.square().clone());
        }
    }
}
