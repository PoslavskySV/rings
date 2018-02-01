package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.WithVariables;
import cc.redberry.rings.poly.MultivariateRing;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static cc.redberry.rings.poly.multivar.GroebnerBasis.canonicalize;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;

/**
 * Ideal represented by its Groebner basis.
 *
 * @since 2.3
 */
public final class Ideal<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements Iterable<Poly>, WithVariables, Serializable {
    /** list of original generators */
    private final List<Poly> originalGenerators;
    /** monomial order used for standard basis */
    private final Comparator<DegreeVector> monomialOrder;
    /** util factory polynomial (ordered by monomialOrder) */
    private final Poly factory;
    /** Groebner basis with respect to {@code monomialOrder} */
    private final Poly[] groebnerBasis;
    /** the whole ring instance (ordered by monomialOrder) */
    private final MultivariateRing<Poly> ring;

    private Ideal(List<Poly> originalGenerators, List<Poly> groebnerBasis) {
        this.originalGenerators = Collections.unmodifiableList(originalGenerators);
        this.factory = groebnerBasis.get(0).createZero();
        this.groebnerBasis = groebnerBasis.toArray(factory.createArray(groebnerBasis.size()));
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
        // todo use FGLM or Groebner walk to switch order
        return monomialOrder == newMonomialOrder
                ? this
                : Ideal(originalGenerators, newMonomialOrder);
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

    Poly mod0(Poly poly) {
        return MultivariateDivision.remainder(setOrdering(poly), groebnerBasis);
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
     * Returns i-th basis generator of this ideal
     */
    public Poly generator(int i) {
        return groebnerBasis[i];
    }

    /**
     * Returns the number of basis generators
     */
    public int nGenerators() {
        return groebnerBasis.length;
    }

    /**
     * Whether this ideal is the whole ring (basis consists of pne constant polynomial)
     */
    public boolean isTrivial() {
        return nGenerators() == 1 && generator(0).isConstant() && !generator(0).isZero();
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
        return nGenerators() == 1 && generator(0).isZero();
    }

    /**
     * Whether this ideal is principal
     */
    public boolean isPrincipal() {
        return nGenerators() == 1;
    }

    /**
     * Whether this ideal contains the specified one
     */
    public boolean contains(Ideal<Term, Poly> oth) {
        return quotient(oth).isTrivial();
    }

    /**
     * Whether this ideal contains the prooduct of two specified ideals
     */
    public boolean containsProduct(Ideal<Term, Poly> a, Ideal<Term, Poly> b) {
        if (a.groebnerBasis.length > b.groebnerBasis.length)
            return containsProduct(b, a);
        return quotient(a).contains(b);
    }

    /**
     * Tests whether specified poly is an element of this ideal
     */
    public boolean contains(Poly poly) {
        return mod0(poly).isZero();
    }

    /**
     * Tests whether {@code poly} belongs to the radical of this
     */
    public boolean isRadicalMember(Poly poly) {
        // adjoin new variable to all generators (convert to F[X][y])
        List<Poly> yGenerators = Arrays.stream(groebnerBasis)
                .map(AMultivariatePolynomial::joinNewVariable)
                .collect(Collectors.toList());

        Poly yPoly = poly.joinNewVariable();
        // add 1 - y*poly
        yGenerators.add(yPoly.createOne().subtract(yPoly.createMonomial(yPoly.nVariables - 1, 1).multiply(yPoly)));
        return Ideal(yGenerators).isTrivial();
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
        l.addAll(Arrays.asList(groebnerBasis));
        l.addAll(Arrays.asList(oth.groebnerBasis));
        return Ideal(l, monomialOrder);
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
        return Ideal(generators, monomialOrder);
    }

    /**
     * Returns the product of this and oth
     */
    public Ideal<Term, Poly> multiply(Poly oth) {
        factory.assertSameCoefficientRingWith(oth);

        if (isTrivial())
            return Ideal(Collections.singletonList(oth), monomialOrder);
        if (oth.isZero())
            return trivial(oth, monomialOrder);
        if (oth.isOne() || this.isEmpty())
            return this;

        return new Ideal<>(canonicalize(
                Arrays.stream(groebnerBasis)
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
            return Ideal(Collections.singletonList(ring.lcm(generator(0), oth.generator(0))), monomialOrder);

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
        tGenerators = GroebnerBasis.GroebnerBasis(tGenerators, blockOrder);
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
        return Ideal(intersection(Ideal(oth)).stream().map(p -> ring.quotient(p, oth)).collect(Collectors.toList()));
    }

    /**
     * Returns the quotient this : oth
     */
    public Ideal<Term, Poly> quotient(Ideal<Term, Poly> oth) {
        if (oth.isEmpty())
            return trivial(factory);
        if (oth.isTrivial())
            return this;
        return oth.stream().map(this::quotient).reduce(trivial(factory), Ideal::intersection);
    }

    Ideal<Term, Poly> insertVariable(int variable) {
        return new Ideal<>(Arrays.stream(groebnerBasis).map(p -> p.insertVariable(variable)).collect(Collectors.toList()));
    }

    @Override
    public Iterator<Poly> iterator() {
        return Arrays.asList(groebnerBasis).iterator();
    }

    /**
     * Stream of generators
     */
    public Stream<Poly> stream() {
        return Arrays.stream(groebnerBasis);
    }

    private void assertSameDomain(Ideal<Term, Poly> oth) {
        factory.assertSameCoefficientRingWith(oth.factory);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Ideal<?, ?> ideal = (Ideal<?, ?>) o;
        // todo: wtd7
        // if (monomialOrder != ideal.monomialOrder)
        //    ideal = ideal.setMonomialOrder(monomialOrder);

        return Arrays.equals(groebnerBasis, ideal.groebnerBasis);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(groebnerBasis);
    }

    @Override
    public String toString(String[] variables) {
        return "[" + stream().map(p -> p.toString(variables)).collect(Collectors.joining(", ")) + "]";
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(factory.nVariables));
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to GREVLEX order will be used.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> Ideal(List<Poly> generators) {
        return Ideal(generators, GREVLEX);
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to GREVLEX order will be used.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> Ideal(Poly... generators) {
        return Ideal(Arrays.asList(generators));
    }

    /**
     * Creates ideal given by a list of generators. Groebner basis with respect to specified {@code monomialOrder} will
     * be used.
     *
     * @param monomialOrder monomial order for unique Groebner basis of the ideal
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> Ideal(List<Poly> generators, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(generators, GroebnerBasis.GroebnerBasis(generators, monomialOrder));
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
        return Ideal(Arrays.stream(generators).
                map(p -> MultivariatePolynomial.parse(p, field, monomialOrder, variables))
                .collect(Collectors.toList()), monomialOrder);
    }
}
