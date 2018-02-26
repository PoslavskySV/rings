package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.WithVariables;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.multivar.GroebnerBasis.*;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static cc.redberry.rings.poly.multivar.GroebnerBasis.GroebnerBasis;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;

/**
 * Ideal represented by its Groebner basis.
 *
 * @since 2.3
 */
public final class Ideal<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements WithVariables, Serializable {
    /** list of original generators */
    private final List<Poly> originalGenerators;
    /** monomial order used for standard basis */
    private final Comparator<DegreeVector> monomialOrder;
    /** util factory polynomial (ordered by monomialOrder) */
    private final Poly factory;
    /** Groebner basis with respect to {@code monomialOrder} */
    private final List<Poly> groebnerBasis;
    /** the whole ring instance (ordered by monomialOrder) */
    private final MultivariateRing<Poly> ring;

    private Ideal(List<Poly> originalGenerators, List<Poly> groebnerBasis) {
        this.originalGenerators = Collections.unmodifiableList(originalGenerators);
        this.factory = groebnerBasis.get(0).createZero();
        this.groebnerBasis = groebnerBasis;
        this.monomialOrder = factory.ordering;
        this.ring = Rings.MultivariateRing(factory);
    }

    private Ideal(List<Poly> groebnerBasis) {
        this(groebnerBasis, groebnerBasis);
    }

    /**
     * The monomial order used for Groebner basis
     */
    public Comparator<DegreeVector> getMonomialOrder() {
        return monomialOrder;
    }

    /**
     * Set the monomial order used for Groebner basis of this ideal
     */
    public Ideal<Term, Poly> setMonomialOrder(Comparator<DegreeVector> newMonomialOrder) {
        if (monomialOrder == newMonomialOrder)
            return this;
        if (isGradedOrder(monomialOrder) || !isGradedOrder(newMonomialOrder))
            return new Ideal<>(originalGenerators, HilbertConvertBasis(groebnerBasis, newMonomialOrder));

        return create(originalGenerators, newMonomialOrder);
    }

    @SuppressWarnings("unchecked")
    private static <Poly extends AMultivariatePolynomial<?, Poly>>
    Poly setOrdering(Poly poly, Comparator<DegreeVector> monomialOrder) {
        return poly.ordering == monomialOrder
                ? poly
                : poly.setOrdering(monomialOrder);
    }

    /** set ordering of poly to monomialOrder */
    private Poly setOrdering(Poly poly) {
        return setOrdering(poly, monomialOrder);
    }

    private Poly mod0(Poly poly) {
        return MultivariateDivision.pseudoRemainder(setOrdering(poly), groebnerBasis);
    }

    /**
     * Reduces {@code poly} modulo this ideal
     */
    public Poly mod(Poly poly) {
        Comparator<DegreeVector> originalOrder = poly.ordering;
        return setOrdering(mod0(poly), originalOrder);
    }

    /**
     * Returns the list of original generators
     */
    public List<Poly> getOriginalGenerators() {
        return originalGenerators;
    }

    /**
     * Groebner basis of this ideal
     */
    public List<Poly> getGroebnerBasis() {
        return Collections.unmodifiableList(groebnerBasis);
    }

    /**
     * Returns the number of elements in Groebner basis
     */
    public int nBasisGenerators() {
        return groebnerBasis.size();
    }

    /**
     * Returns i-th element of Groebner basis
     */
    public Poly getBasisGenerator(int i) {
        return groebnerBasis.get(i);
    }

    /**
     * Whether this ideal is the whole ring (basis consists of pne constant polynomial)
     */
    public boolean isTrivial() {
        return nBasisGenerators() == 1 && getBasisGenerator(0).isConstant() && !getBasisGenerator(0).isZero();
    }

    /**
     * Whether this is a proper ideal
     */
    public boolean isProper() {
        return !isTrivial();
    }

    /**
     * Whether this ideal is empty
     */
    public boolean isEmpty() {
        return nBasisGenerators() == 1 && getBasisGenerator(0).isZero();
    }

    /**
     * Whether this ideal is principal
     */
    public boolean isPrincipal() {
        return nBasisGenerators() == 1;
    }

    /**
     * Whether this ideal is homogeneous
     */
    public boolean isHomogeneous() {
        return isHomogeneousIdeal(groebnerBasis);
    }

    /**
     * Whether this ideal is monomial
     */
    public boolean isMonomial() {
        return isMonomialIdeal(groebnerBasis);
    }

    /**
     * Tests whether specified poly is an element of this ideal
     */
    public boolean contains(Poly poly) {
        return mod0(poly).isZero();
    }

    /**
     * Whether this ideal contains the specified one
     */
    public boolean contains(Ideal<Term, Poly> oth) {
        return quotient(oth).isTrivial();
    }

    // lazy Hilbert-Poincare series
    private HilbertSeries hilbertSeries = null;

    /** Hilbert-Poincare series of this ideal */
    public synchronized HilbertSeries getHilbertSeries() {
        if (hilbertSeries == null) {
            if (isHomogeneous() || isGradedOrder(monomialOrder))
                hilbertSeries = HilbertSeriesOfLeadingTermsSet(groebnerBasis);
            else
                // use original generators to construct basis when current ordering is "hard"
                hilbertSeries = HilbertSeriesOfLeadingTermsSet(
                        GroebnerBasisWithOptimizedGradedOrder(originalGenerators));
        }
        return hilbertSeries;
    }

    /** Returns the affine dimension of this ideal */
    public int dimension() {
        return getHilbertSeries().dimension();
    }

    /** Returns the affine degree of this ideal */
    public int degree() {
        return getHilbertSeries().degree();
    }

    /**
     * Whether this ideal contains the product of two specified ideals
     */
    public boolean containsProduct(Ideal<Term, Poly> a, Ideal<Term, Poly> b) {
        if (a.nBasisGenerators() > b.nBasisGenerators())
            return containsProduct(b, a);
        return quotient(a).contains(b);
    }

    /**
     * Tests whether {@code poly} belongs to the radical of this
     */
    public boolean isRadicalMember(Poly poly) {
        // adjoin new variable to all generators (convert to F[X][y])
        List<Poly> yGenerators = groebnerBasis.stream()
                .map(AMultivariatePolynomial::joinNewVariable)
                .collect(Collectors.toList());

        Poly yPoly = poly.joinNewVariable();
        // add 1 - y*poly
        yGenerators.add(yPoly.createOne().subtract(yPoly.createMonomial(yPoly.nVariables - 1, 1).multiply(yPoly)));
        return create(yGenerators).isTrivial();
    }

    /**
     * Returns the union of this and oth
     */
    public Ideal<Term, Poly> union(Poly oth) {
        factory.assertSameCoefficientRingWith(oth);
        if (oth.isZero())
            return this;
        if (oth.isOne())
            return trivial(factory);

        List<Poly> l = new ArrayList<>(groebnerBasis);
        l.add(oth);
        return create(l, monomialOrder);
    }

    /**
     * Returns the union of this and oth
     */
    public Ideal<Term, Poly> union(Ideal<Term, Poly> oth) {
        assertSameDomain(oth);
        if (isEmpty() || oth.isTrivial())
            return oth;
        if (oth.isEmpty() || isTrivial())
            return this;

        List<Poly> l = new ArrayList<>();
        l.addAll(groebnerBasis);
        l.addAll(oth.groebnerBasis);
        return create(l, monomialOrder);
    }

    /**
     * Returns the product of this and oth
     */
    public Ideal<Term, Poly> multiply(Ideal<Term, Poly> oth) {
        assertSameDomain(oth);
        if (isTrivial() || oth.isEmpty())
            return oth;
        if (oth.isTrivial() || this.isEmpty())
            return this;

        List<Poly> generators = new ArrayList<>();
        for (Poly a : groebnerBasis)
            for (Poly b : oth.groebnerBasis)
                generators.add(a.clone().multiply(b));
        return create(generators, monomialOrder);
    }

    /**
     * Returns squared ideal
     */
    public Ideal<Term, Poly> square() {
        return multiply(this);
    }

    /**
     * Returns this in a power of exponent
     */
    public Ideal<Term, Poly> pow(int exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        if (exponent == 1)
            return this;

        Ideal<Term, Poly> result = trivial(factory);
        Ideal<Term, Poly> k2p = this;
        for (; ; ) {
            if ((exponent & 1) != 0)
                result = result.multiply(k2p);
            exponent >>= 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p);
        }
    }

    /**
     * Returns the product of this and oth
     */
    public Ideal<Term, Poly> multiply(Poly oth) {
        factory.assertSameCoefficientRingWith(oth);

        if (isTrivial())
            return create(Collections.singletonList(oth), monomialOrder);
        if (oth.isZero())
            return trivial(oth, monomialOrder);
        if (oth.isOne() || this.isEmpty())
            return this;

        return new Ideal<>(canonicalize(
                groebnerBasis.stream()
                        .map(p -> p.clone().multiply(oth))
                        .collect(Collectors.toList())));
    }

    /**
     * Returns the intersection of this and oth
     */
    @SuppressWarnings("unchecked")
    public Ideal<Term, Poly> intersection(Ideal<Term, Poly> oth) {
        assertSameDomain(oth);

        if (isTrivial() || oth.isEmpty())
            return oth;
        if (oth.isTrivial() || this.isEmpty())
            return this;
        if (isPrincipal() && oth.isPrincipal())
            // intersection of principal ideals is easy
            return create(Collections.singletonList(ring.lcm(getBasisGenerator(0), oth.getBasisGenerator(0))), monomialOrder);

        // we compute (t * I + (1 - t) * J) ∩ R[X]
        Poly t = factory.insertVariable(0).createMonomial(0, 1);
        List<Poly> tGenerators = new ArrayList<>();
        for (Poly gI : this.groebnerBasis)
            tGenerators.add(gI.insertVariable(0).multiply(t));
        Poly omt = t.clone().negate().increment(); // 1 - t
        for (Poly gJ : oth.groebnerBasis)
            tGenerators.add(gJ.insertVariable(0).multiply(omt));

        Comparator<DegreeVector> blockOrder = MonomialOrder.product(
                MonomialOrder.LEX, 1,
                monomialOrder, factory.nVariables);

        // elimination
        tGenerators = GroebnerBasis(tGenerators, blockOrder);
        List<Poly> result = tGenerators.stream()
                .filter(p -> p.degree(0) == 0)
                .map(p -> p.dropVariable(0))
                .map(p -> p.setOrdering(monomialOrder)) // <- restore order!
                .collect(Collectors.toList());
        canonicalize(result);
        return new Ideal<>(result);
    }

    /**
     * Returns the quotient this : oth
     */
    @SuppressWarnings("unchecked")
    public Ideal<Term, Poly> quotient(Poly oth) {
        if (oth.isZero())
            return trivial(factory);
        if (oth.isConstant())
            return this;
        return create(intersection(create(oth)).groebnerBasis.stream().map(p -> ring.quotient(p, oth)).collect(Collectors.toList()));
    }

    /**
     * Returns the quotient this : oth
     */
    public Ideal<Term, Poly> quotient(Ideal<Term, Poly> oth) {
        if (oth.isEmpty())
            return trivial(factory);
        if (oth.isTrivial())
            return this;
        return oth.groebnerBasis.stream().map(this::quotient).reduce(trivial(factory), Ideal::intersection);
    }

    Ideal<Term, Poly> insertVariable(int variable) {
        return new Ideal<>(groebnerBasis.stream().map(p -> p.insertVariable(variable)).collect(Collectors.toList()));
    }

    private void assertSameDomain(Ideal<Term, Poly> oth) {
        factory.assertSameCoefficientRingWith(oth.factory);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Ideal<?, ?> ideal = (Ideal<?, ?>) o;
        return monomialOrder.equals(ideal.monomialOrder)
                && groebnerBasis.equals(ideal.groebnerBasis);
    }

    @Override
    public int hashCode() {
        return groebnerBasis.hashCode();
    }

    @Override
    public String toString(String[] variables) {
        return "<" + groebnerBasis.stream().map(p -> p.toString(variables)).collect(Collectors.joining(", ")) + ">";
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(factory.nVariables));
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to GREVLEX order will be used.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> create(List<Poly> generators) {
        return create(generators, GREVLEX);
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to GREVLEX order will be used.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> create(Poly... generators) {
        return create(Arrays.asList(generators));
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to specified {@code monomialOrder} will
     * be used.
     *
     * @param monomialOrder monomial order for unique Groebner basis of the ideal
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> create(List<Poly> generators, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(generators, GroebnerBasis(generators, monomialOrder));
    }

    /**
     * Creates trivial ideal (ideal = ring)
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> trivial(Poly factory) {
        return trivial(factory, GREVLEX);
    }

    /**
     * Creates trivial ideal (ideal = ring)
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> trivial(Poly factory, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(Collections.singletonList(factory.createOne().setOrdering(monomialOrder)));
    }

    /**
     * Creates empty ideal
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> empty(Poly factory) {
        return empty(factory, GREVLEX);
    }

    /**
     * Creates empty ideal
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> empty(Poly factory, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(Collections.singletonList(factory.createZero().setOrdering(monomialOrder)));
    }

    /**
     * Shortcut for parse
     */
    public static <E> Ideal<Monomial<E>, MultivariatePolynomial<E>>
    parse(String[] generators, Ring<E> field, String[] variables) {
        return parse(generators, field, GREVLEX, variables);
    }

    /**
     * Shortcut for parse
     */
    public static <E> Ideal<Monomial<E>, MultivariatePolynomial<E>>
    parse(String[] generators, Ring<E> field, Comparator<DegreeVector> monomialOrder, String[] variables) {
        return create(Arrays.stream(generators).
                map(p -> MultivariatePolynomial.parse(p, field, monomialOrder, variables))
                .collect(Collectors.toList()), monomialOrder);
    }
}
