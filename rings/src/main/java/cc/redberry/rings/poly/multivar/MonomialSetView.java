package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;

import java.util.Collection;
import java.util.Iterator;

/**
 * Parent auxiliary interface for multivariate polynomials.
 *
 * @since 2.3
 */
interface MonomialSetView<Term extends AMonomial<Term>> extends Iterable<Term> {
    Iterator<Term> ascendingIterator();

    Iterator<Term> descendingIterator();

    @Override
    default Iterator<Term> iterator() { return ascendingIterator(); }

    default Term first() { return ascendingIterator().next(); }

    default Term last() { return descendingIterator().next(); }

    default Term lt() { return descendingIterator().next(); }

    /**
     * Returns an array of degrees of all variables, so that is i-th element of the result is the polynomial degree
     * (univariate) with respect to i-th variable
     *
     * @return array of degrees
     */
    int[] degrees();

    /**
     * Returns the sum of {@link #degrees()}
     *
     * @return sum of {@link #degrees()}
     */
    default int degreeSum() { return ArraysUtil.sum(degrees());}

    /**
     * Collection view of all terms
     */
    Collection<Term> collection();
}
