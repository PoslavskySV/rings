package cc.redberry.rings.poly;

import cc.redberry.rings.Factors;
import cc.redberry.rings.ToStringSupport;
import cc.redberry.rings.WithVariables;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.PolynomialRing;
import static cc.redberry.rings.poly.PolynomialMethods.polyPow;

/**
 * {@inheritDoc}
 *
 * @since 1.0
 */
public final class FactorDecomposition<Poly extends IPolynomial<Poly>>
        extends Factors<Poly>
        implements WithVariables, java.io.Serializable {
    private static final long serialVersionUID = 1L;

    private FactorDecomposition(Poly unit, List<Poly> factors, TIntArrayList exponents) {
        super(PolynomialRing(unit), unit, factors, exponents);
    }

    private FactorDecomposition(Factors<Poly> factors) {
        super(factors.ring, factors.unit, factors.factors, factors.exponents);
    }

    @Override
    public boolean isUnit(Poly element) {
        return element.isConstant();
    }

    @Override
    public FactorDecomposition<Poly> setUnit(Poly unit) {
        super.setUnit(unit);
        return this;
    }

    @Override
    public FactorDecomposition<Poly> addUnit(Poly unit) {
        super.addUnit(unit);
        return this;
    }

    @Override
    public FactorDecomposition<Poly> addFactor(Poly factor, int exponent) {
        super.addFactor(factor, exponent);
        return this;
    }

    @Override
    public FactorDecomposition<Poly> addAll(Factors<Poly> other) {
        super.addAll(other);
        return this;
    }

    @Override
    public FactorDecomposition<Poly> canonical() {
        if (factors.size() == 0)
            return this;
        reduceUnitContent();
        Poly[] fTmp = factors.toArray(factors.get(0).createArray(factors.size()));
        int[] eTmp = exponents.toArray();
        for (int i = fTmp.length - 1; i >= 0; --i) {
            Poly poly = fTmp[i];
            if (poly.isMonomial() && eTmp[i] != 1) {
                poly = PolynomialMethods.polyPow(poly, eTmp[i], false);
                assert poly.isMonomial();
            }
            if (poly.signumOfLC() < 0) {
                poly.negate();
                if (eTmp[i] % 2 == 1)
                    unit.negate();
            }
        }

        ArraysUtil.quickSort(fTmp, eTmp);
        for (int i = 0; i < fTmp.length; i++) {
            factors.set(i, fTmp[i]);
            exponents.set(i, eTmp[i]);
        }
        return this;
    }

    /**
     * Calculates the signum of the polynomial constituted by this decomposition
     *
     * @return the signum of the polynomial constituted by this decomposition
     */
    public int signum() {
        int signum = unit.signumOfLC();
        for (int i = 0; i < factors.size(); i++)
            signum *= exponents.get(i) % 2 == 0 ? 1 : factors.get(i).signumOfLC();
        return signum;
    }

    /**
     * Makes each factor monic (moving leading coefficients to the {@link #unit})
     */
    public FactorDecomposition<Poly> monic() {
        for (int i = 0; i < factors.size(); i++) {
            Poly factor = factors.get(i);
            addUnit(polyPow(factor.lcAsPoly(), exponents.get(i), false));
            factor = factor.monic();
            assert factor != null;
        }
        return this;
    }

    /**
     * Makes each factor primitive (moving contents to the {@link #unit})
     */
    public FactorDecomposition<Poly> primitive() {
        for (int i = 0; i < factors.size(); i++) {
            Poly factor = factors.get(i);
            Poly content = factor.contentAsPoly();
            addUnit(polyPow(content, exponents.get(i), false));
            factor = factor.divideByLC(content);
            assert factor != null;
            if (factor.signumOfLC() < 0) {
                factor.negate();
                if (exponents.get(i) % 2 == 1)
                    unit.negate();
            }
        }
        return this;
    }

    public <OthPoly extends IPolynomial<OthPoly>> FactorDecomposition<OthPoly> map(Function<Poly, OthPoly> mapper) {
        return of(mapper.apply(unit), factors.stream().map(mapper).collect(Collectors.toList()), exponents);
    }

    /**
     * Calls {@link #monic()} if the coefficient ring is field and {@link #primitive()} otherwise
     */
    public FactorDecomposition<Poly> reduceUnitContent() {
        return unit.isOverField() ? monic() : primitive();
    }

    @Override
    public String toString(String[] variables) {
        return toString(ToStringSupport.withVariables(variables));
    }

    @Override
    public FactorDecomposition<Poly> clone() {
        return new FactorDecomposition<>(unit.clone(), factors.stream().map(Poly::clone).collect(Collectors.toList()), new TIntArrayList(exponents));
    }

    /** Unit factorization */
    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly> unit(Poly unit) {
        if (!unit.isConstant())
            throw new IllegalArgumentException();
        return empty(unit).addUnit(unit);
    }

    /** Empty factorization */
    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly> empty(Poly factory) {
        return new FactorDecomposition<>(factory.createOne(), new ArrayList<>(), new TIntArrayList());
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param unit      the unit coefficient
     * @param factors   the factors
     * @param exponents the exponents
     */
    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly>
    of(Poly unit, List<Poly> factors, TIntArrayList exponents) {
        if (factors.size() != exponents.size())
            throw new IllegalArgumentException();
        FactorDecomposition<Poly> r = empty(unit).addUnit(unit);
        for (int i = 0; i < factors.size(); i++)
            r.addFactor(factors.get(i), exponents.get(i));
        return r;
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param factors factors
     */
    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly>
    of(Poly... factors) {
        if (factors.length == 0)
            throw new IllegalArgumentException();
        return of(Arrays.asList(factors));
    }

    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly>
    of(Poly a) {
        Poly[] array = a.createArray(1);
        array[0] = a;
        return of(array);
    }

    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly>
    of(Poly a, Poly b) {
        return of(a.createArray(a, b));
    }

    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly>
    of(Poly a, Poly b, Poly c) {
        return of(a.createArray(a, b, c));
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param factors factors
     */
    public static <Poly extends IPolynomial<Poly>> FactorDecomposition<Poly> of(Collection<Poly> factors) {
        TObjectIntHashMap<Poly> map = new TObjectIntHashMap<>();
        for (Poly e : factors)
            map.adjustOrPutValue(e, 1, 1);
        List<Poly> l = new ArrayList<>();
        TIntArrayList e = new TIntArrayList();
        map.forEachEntry((a, b) -> {
            l.add(a);
            e.add(b);
            return true;
        });
        return of(factors.iterator().next().createOne(), l, e);
    }
}
