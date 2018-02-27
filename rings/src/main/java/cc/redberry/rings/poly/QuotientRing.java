package cc.redberry.rings.poly;

import cc.redberry.rings.AQuotientRing;
import cc.redberry.rings.WithVariables;
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
        return ideal.mod(el);
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
    public Poly parse(String string) {
        return parse(string, WithVariables.defaultVars(nVariables()));
    }

    @Override
    public Poly parse(String string, String[] variables) {
        return valueOf(factory.parsePoly(string, variables));
    }

    @Override
    public Poly variable(int variable) {
        return factory.createMonomial(variable, 1);
    }

    @Override
    public String toString(String[] variables) {
        return toString(factory.coefficientRingToString(), variables);
    }

    @Override
    public String toString(String coefficientDomain, String[] variables) {
        return ((IPolynomialRing) baseRing).toString(coefficientDomain, variables) + "/" + ideal.toString(variables) + "";
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(nVariables()));
    }
}
