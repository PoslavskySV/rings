package cc.redberry.rings.poly.multivar;

import java.util.*;

/**
 * Sorted set of monomials -- basic underlying data structure of multivariate polynomials.
 *
 * @since 1.0
 */
public final class MonomialSet<Term extends AMonomial<Term>>
        extends TreeMap<DegreeVector, Term>
        implements MonomialSetView<Term>, Iterable<Term>, Cloneable {
    private static final long serialVersionUID = 1L;

    public MonomialSet(Comparator<? super DegreeVector> comparator) {
        super(comparator);
    }

    /**
     * Constructs a new monomial set containing the same mappings and using the same ordering as the specified sorted
     * map.  This method runs in linear time.
     *
     * @param m the sorted map whose mappings are to be placed in this monomial set, and whose comparator is to be used
     *          to sort this map
     */
    public MonomialSet(SortedMap<DegreeVector, ? extends Term> m) {
        super(m);
    }

    @Override
    public Iterator<Term> iterator() {
        return values().iterator();
    }

    /**
     * Add monomial to this set
     *
     * @param term monomial
     * @return this
     */
    public Term add(Term term) {return put(term, term);}

    /** First monomial in this set */
    @Override
    public Term first() {return firstEntry().getValue();}

    /** Last monomial in this set */
    @Override
    public Term last() {return lastEntry().getValue();}

    @Override
    public Iterator<Term> ascendingIterator() {
        return values().iterator();
    }

    @Override
    public Iterator<Term> descendingIterator() {
        return descendingMap().values().iterator();
    }

    @Override
    public Collection<Term> collection() {
        return values();
    }

    @Override
    public int[] degrees() {
        throw new UnsupportedOperationException();
    }

    @Override
    @SuppressWarnings("unchecked")
    public MonomialSet<Term> clone() {
        return (MonomialSet<Term>) super.clone();
    }

    @Override
    public int hashCode() {
        int h = 0;
        Iterator<Map.Entry<DegreeVector, Term>> i = entrySet().iterator();
        while (i.hasNext())
            h += i.next().getValue().hashCode();
        return h;
    }

    public int skeletonHashCode() {
        int h = 0;
        Iterator<Map.Entry<DegreeVector, Term>> i = entrySet().iterator();
        while (i.hasNext())
            h += i.next().getKey().dv().hashCode();
        return h;
    }
}
