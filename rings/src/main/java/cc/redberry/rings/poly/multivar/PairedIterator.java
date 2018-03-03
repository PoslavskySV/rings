package cc.redberry.rings.poly.multivar;

import java.util.Comparator;
import java.util.Iterator;

/**
 * Iterator over a pair of polynomials
 */
public final class PairedIterator<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
    final Poly factory;
    final Term zeroTerm;
    final Comparator<DegreeVector> ordering;
    final Iterator<Term> aIterator, bIterator;

    public PairedIterator(Poly a, Poly b) {
        this.factory = a;
        this.zeroTerm = factory.monomialAlgebra.getZeroTerm(factory.nVariables);
        this.ordering = a.ordering;
        this.aIterator = a.iterator();
        this.bIterator = b.iterator();
    }

    public boolean hasNext() {
        return aIterator.hasNext() || bIterator.hasNext() || aTermCached != null || bTermCached != null;
    }

    public Term aTerm = null, bTerm = null;
    private Term aTermCached = null, bTermCached = null;

    public void advance() {
        if (aTermCached != null) {
            aTerm = aTermCached;
            aTermCached = null;
        } else
            aTerm = aIterator.hasNext() ? aIterator.next() : zeroTerm;
        if (bTermCached != null) {
            bTerm = bTermCached;
            bTermCached = null;
        } else
            bTerm = bIterator.hasNext() ? bIterator.next() : zeroTerm;

        int c = ordering.compare(aTerm, bTerm);
        if (c < 0 && aTerm != zeroTerm) {
            bTermCached = bTerm;
            bTerm = zeroTerm;
        } else if (c > 0 && bTerm != zeroTerm) {
            aTermCached = aTerm;
            aTerm = zeroTerm;
        }
    }
}
