package cc.r2.core.poly.multivar2;

import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.set.hash.TIntHashSet;

import java.util.Comparator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
abstract class AbstractMultivariatepolynomial<
        Monomial extends DegreeVector<Monomial>,
        ThisType extends AbstractMultivariatepolynomial<Monomial, ThisType>>
        implements IGeneralPolynomial<ThisType> {
    final int nVariables;
    final Comparator<DegreeVector> ordering;
    final MonomialsSet<Monomial> data;

    AbstractMultivariatepolynomial(int nVariables, Comparator<DegreeVector> ordering, MonomialsSet<Monomial> data) {
        this.nVariables = nVariables;
        this.ordering = ordering;
        this.data = data;
    }

    abstract ThisType create(MonomialsSet<Monomial> set);

    abstract ThisType add(Monomial monomial);

    /** release caches */
    void release() {}

    /**
     * Copies this with the new ordering {@code newOrdering}
     *
     * @param newOrdering the new ordering
     * @return a copy of this with a new ordering
     */
    public ThisType setOrdering(Comparator<DegreeVector> newOrdering) {
        MonomialsSet<Monomial> newData = new MonomialsSet<>(newOrdering);
        newData.putAll(data);
        return create(newData);
    }

    /**
     * Returns the number of terms in this polynomial
     *
     * @return number of terms
     */
    public int size() {return data.size();}

    @Override
    public boolean isZero() {
        return data.size() == 0;
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && data.first().isZeroVector());
    }

    @Override
    public boolean isMonomial() {
        return size() <= 1;
    }

    /**
     * Returns whether this poly is effectively univariate
     *
     * @return whether this poly is effectively univariate
     */
    public boolean isEffectiveUnivariate() {
        int[] degrees = degrees();
        boolean b = false;
        for (int i = 0; i < nVariables; i++) {
            if (degrees[i] != 0) {
                if (b)
                    return false;
                else
                    b = true;
            }
        }
        return true;
    }

    /**
     * Returns the total degree of this polynomial, that is the maximal degree among all terms
     *
     * @return the total degree of this polynomial, that is the maximal degree among all terms
     */
    @Override
    public int degree() {
        int max = 0;
        for (Monomial db : data)
            max = Math.max(max, db.totalDegree);
        return max;
    }

    /**
     * Returns the degree of this polynomial with respect to i-th variable
     *
     * @return the degree of this polynomial with respect to i-th variable
     */
    public int degree(int i) {
        int max = 0;
        for (Monomial db : data)
            max = Math.max(max, db.exponents[i]);
        return max;
    }

    /**
     * Returns an array of degrees of all variables, so that is i-th element of the result is the polynomial degree
     * with respect to i-th variable
     *
     * @return array of degrees
     */
    public int[] degrees() {
        int[] degrees = new int[nVariables];
        for (Monomial db : data)
            for (int i = 0; i < nVariables; i++)
                if (db.exponents[i] > degrees[i])
                    degrees[i] = db.exponents[i];
        return degrees;
    }

    /**
     * Returns the degrees in which {@code variable} occurs in this polynomial
     *
     * @return the degrees in which {@code variable} occurs in this polynomial
     */
    public int[] degrees(int variable) {
        TIntHashSet degrees = new TIntHashSet();
        for (Monomial db : data)
            degrees.add(db.exponents[variable]);
        return degrees.toArray();
    }

    /**
     * Returns the sum of {@link #degrees()}
     *
     * @return sum of {@link #degrees()}
     */
    public int degreeSum() {
        int r = 0;
        for (int d : degrees())
            r += d;
        return r;
    }

    /**
     * Returns the number of really used variables
     *
     * @return the number of presenting variables
     */
    public int usedVariables() {
        int[] degrees = degrees();
        int r = 0;
        for (int d : degrees)
            if (d != 0)
                ++r;
        return r;
    }

    /** check whether number of variables is the same */
    void ensureCompatible(ThisType oth) {
        if (nVariables != oth.nVariables)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /** check whether number of variables is the same */
    void ensureCompatible(Monomial oth) {
        if (nVariables != oth.exponents.length)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    @SuppressWarnings("unchecked")
    private ThisType self = (ThisType) this;

    @Override
    public ThisType toZero() {
        data.clear();
        release();
        return self;
    }

    @Override
    public ThisType set(ThisType oth) {
        ensureCompatible(oth);
        data.clear();
        data.putAll(oth.data);
        release();
        return self;
    }

    @Override
    public ThisType createZero() {
        return create(new MonomialsSet<>(ordering));
    }

    /**
     * Returns a copy of this with {@code nVariables + 1}
     *
     * @return a copy of this with with one additional (last) variable added
     */
    public ThisType joinNewVariable() {
        MonomialsSet<Monomial> newData = new MonomialsSet<>(ordering);
        for (Monomial term : data)
            newData.add(term.joinNewVariable());
        return create(newData);
    }

    /**
     * Returns a coefficient before {@code variable^exponent} as a multivariate polynomial
     *
     * @param variable the variable
     * @param exponent the exponent
     * @return coefficient before {@code variable^exponent} as a multivariate polynomial
     */
    public ThisType coefficientOf(int variable, int exponent) {
        ThisType result = createZero();
        for (Monomial e : data) {
            if (e.exponents[variable] != exponent)
                continue;
            result.add(e.setZero(variable));
        }
        return result;
    }

    /**
     * Renames variable {@code i} to {@code j} and {@code j} to {@code i}
     *
     * @param i the first variable
     * @param j the second variable
     * @return polynomial with variable {@code i} renamed to {@code j} and {@code j} renamed to {@code i}
     */
    public ThisType swapVariables(int i, int j) {
        int[] newVariables = ArraysUtil.sequence(nVariables);
        newVariables[i] = j;
        newVariables[j] = i;
        return renameVariables(newVariables);
    }

    /**
     * Rename variables from [0,1,...N] to [newVariables[0], newVariables[1], ..., newVariables[N]]
     *
     * @param newVariables the new variables
     * @return renamed polynomial
     */
    public ThisType renameVariables(int[] newVariables) {
        // NOTE: always return a copy of poly, even if order of variables is unchanged
        MonomialsSet<Monomial> newData = new MonomialsSet<>(ordering);
        for (Monomial e : this.data)
            newData.add(e.setDegreeVector(map(e.exponents, newVariables), e.totalDegree));
        return create(newData);
    }

    @Override
    public abstract  ThisType clone();

    private static int[] map(int[] degrees, int[] mapping) {
        int[] newDegrees = new int[degrees.length];
        for (int i = 0; i < degrees.length; i++)
            newDegrees[i] = degrees[mapping[i]];
        return newDegrees;
    }
}
