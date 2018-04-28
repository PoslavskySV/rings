package cc.redberry.rings.poly;

import cc.redberry.rings.AQuotientRing;
import cc.redberry.rings.io.IStringifier;
import cc.redberry.rings.poly.multivar.AMonomial;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.Ideal;

/**
 * Multivariate quotient ring
 */
public class QuotientRing<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        extends AQuotientRing<Poly> implements IPolynomialRing<Poly> {
    /** the ideal */
    public final Ideal<Term, Poly> ideal;
    private final Poly factory;

    public QuotientRing(IPolynomialRing<Poly> baseRing, Ideal<Term, Poly> ideal) {
        // fixme polys over Z ????
        super(baseRing);
        this.ideal = ideal;
        this.factory = ideal.getBasisGenerator(0).createZero();
    }

    @Override
    public Poly mod(Poly el) {
        return ideal.normalForm(el);
    }

    @Override
    public int nVariables() {
        return factory.nVariables;
    }

    @Override
    public Poly factory() {
        return factory;
    }

    @Override
    public Poly variable(int variable) {
        return factory.createMonomial(variable, 1);
    }


    @Override
    public String toString(IStringifier<Poly> stringifier) {
        return baseRing.toString(stringifier) + "/<" + ideal.toString(stringifier) + ">";
    }

    @Override
    public String toString() {
        return toString(IStringifier.dummy());
    }
}
