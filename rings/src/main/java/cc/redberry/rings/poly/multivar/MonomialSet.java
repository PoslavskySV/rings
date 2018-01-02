package cc.redberry.rings.poly.multivar;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Sorted set of monomials -- basic underlying data structure of multivariate polynomials.
 *
 * @since 1.0
 */
public final class MonomialSet<Term extends AMonomial<Term>>
        extends TreeMap<DegreeVector, Term>
        implements Iterable<Term>, Cloneable {
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
    public Term first() {return firstEntry().getValue();}

    /** Last monomial in this set */
    public Term last() {return lastEntry().getValue();}

    @Override
    @SuppressWarnings("unchecked")
    public MonomialSet<Term> clone() {return (MonomialSet<Term>) super.clone();}
}
