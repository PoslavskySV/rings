package cc.r2.core.poly.multivar;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Sorted set of monomials -- basic underlying data structure of multivariate polynomials.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MonomialSet<Monomial extends DegreeVector>
        extends TreeMap<DegreeVector, Monomial>
        implements Iterable<Monomial>, Cloneable {
    private static final long serialVersionUID = 1L;

    public MonomialSet(Comparator<? super DegreeVector> comparator) {
        super(comparator);
    }

    /**
     * Constructs a new monomial set containing the same mappings and
     * using the same ordering as the specified sorted map.  This
     * method runs in linear time.
     *
     * @param m the sorted map whose mappings are to be placed in this monomial set,
     *          and whose comparator is to be used to sort this map
     */
    public MonomialSet(SortedMap<DegreeVector, ? extends Monomial> m) {
        super(m);
    }

    @Override
    public Iterator<Monomial> iterator() {
        return values().iterator();
    }

    /**
     * Add monomial to this set
     *
     * @param term monomial
     * @return this
     */
    public Monomial add(Monomial term) {return put(term, term);}

    /** First monomial in this set */
    public Monomial first() {return firstEntry().getValue();}

    /** Last monomial in this set */
    public Monomial last() {return lastEntry().getValue();}

    @Override
    @SuppressWarnings("unchecked")
    public MonomialSet<Monomial> clone() {return (MonomialSet<Monomial>) super.clone();}
}
