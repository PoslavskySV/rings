package cc.redberry.rings.poly.multivar;

import java.util.Comparator;
import java.util.Iterator;

/**
 * Iterator over a pair of polynomials
 */
public final class PairedIterator<
        Term1 extends AMonomial<Term1>,
        Poly1 extends AMultivariatePolynomial<Term1, Poly1>,
        Term2 extends AMonomial<Term2>,
        Poly2 extends AMultivariatePolynomial<Term2, Poly2>> {
    final Term1 zeroTerm1;
    final Term2 zeroTerm2;
    final Comparator<DegreeVector> ordering;
    final Iterator<Term1> aIterator;
    final Iterator<Term2> bIterator;

    public PairedIterator(Poly1 a, Poly2 b) {
        this.zeroTerm1 = a.monomialAlgebra.getZeroTerm(a.nVariables);
        this.zeroTerm2 = b.monomialAlgebra.getZeroTerm(b.nVariables);
        this.ordering = a.ordering;
        this.aIterator = a.iterator();
        this.bIterator = b.iterator();
    }

    public boolean hasNext() {
        return aIterator.hasNext() || bIterator.hasNext() || aTermCached != null || bTermCached != null;
    }

    public Term1 aTerm = null;
    public Term2 bTerm = null;
    private Term1 aTermCached = null;
    private Term2 bTermCached = null;

    public void advance() {
        if (aTermCached != null) {
            aTerm = aTermCached;
            aTermCached = null;
        } else
            aTerm = aIterator.hasNext() ? aIterator.next() : zeroTerm1;
        if (bTermCached != null) {
            bTerm = bTermCached;
            bTermCached = null;
        } else
            bTerm = bIterator.hasNext() ? bIterator.next() : zeroTerm2;

        int c = ordering.compare(aTerm, bTerm);
        if (c < 0 && aTerm != zeroTerm1) {
            bTermCached = bTerm;
            bTerm = zeroTerm2;
        } else if (c > 0 && bTerm != zeroTerm2) {
            aTermCached = aTerm;
            aTerm = zeroTerm1;
        }
    }
}
