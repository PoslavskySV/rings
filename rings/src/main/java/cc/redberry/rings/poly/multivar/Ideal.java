package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.WithVariables;
import cc.redberry.rings.poly.MultivariateRing;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Ideal represented by its Groebner basis.
 *
 * @since 2.3
 */
public final class Ideal<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements Iterable<Poly>, WithVariables, Serializable {
    /** util factory polynomial */
    private final Poly factory;
    /** monomial order */
    private final Comparator<DegreeVector> monomialOrder;
    /** Groebner basis */
    private final Poly[] groebnerBasis;
    /** the whole ring instance */
    private final MultivariateRing<Poly> ring;

    private Ideal(Collection<Poly> groebnerBasis) {
        this.factory = groebnerBasis.iterator().next().createZero();
        this.monomialOrder = factory.ordering;
        this.groebnerBasis = groebnerBasis.toArray(factory.createArray(groebnerBasis.size()));
        this.ring = Rings.MultivariateRing(factory);
    }

    /**
     * The monomial order used
     */
    public Comparator<DegreeVector> getMonomialOrder() {
        return monomialOrder;
    }

    /**
     * Reduces {@code poly} modulo this ideal
     */
    public Poly mod(Poly poly) {
        return MultivariateDivision.remainder(poly, groebnerBasis);
    }

    /**
     * Returns i-th generator of this ideal
     */
    public Poly generator(int i) {
        return groebnerBasis[i];
    }

    /**
     * Returns the number of generators
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
     * Tests whether specified poly is an element of this ideal
     */
    public boolean isMember(Poly poly) {
        return mod(poly).isZero();
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
        return Ideal(yGenerators, MonomialOrder.GREVLEX).isTrivial();
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
    public Ideal<Term, Poly> product(Ideal<Term, Poly> oth) {
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
    public Ideal<Term, Poly> product(Poly oth) {
        factory.assertSameCoefficientRingWith(oth);

        Poly othReduced = mod(oth);
        if (isTrivial())
            return Ideal(Collections.singletonList(othReduced), monomialOrder);
        if (othReduced.isZero())
            return trivial(oth, monomialOrder);
        if (othReduced.isOne() || this.isEmpty())
            return this;

        return new Ideal<>(GroebnerBasis.canonicalize(
                Arrays.stream(groebnerBasis)
                        .map(p -> p.clone().multiply(othReduced))
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
                MonomialOrder.GREVLEX, factory.nVariables);

        // elimination
        tGenerators = GroebnerBasis.GroebnerBasis(tGenerators, blockOrder);
        List<Poly> result = tGenerators.stream().filter(p -> p.degree(0) == 0).map(p -> p.dropVariable(0)).collect(Collectors.toList());
        GroebnerBasis.canonicalize(result);
        return new Ideal<>(result);
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
     * Creates ideal given by a list of generators computing corresponding unique Groebner basis
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> Ideal(List<Poly> generators, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(GroebnerBasis.GroebnerBasis(generators, monomialOrder));
    }

    /**
     * Creates trivial ideal (ideal = ring)
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> trivial(Poly factory, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(Collections.singleton(factory.createOne().setOrdering(monomialOrder)));
    }

    /**
     * Creates empty ideal
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Ideal<Term, Poly> empty(Poly factory, Comparator<DegreeVector> monomialOrder) {
        return new Ideal<>(Collections.singleton(factory.createZero().setOrdering(monomialOrder)));
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
