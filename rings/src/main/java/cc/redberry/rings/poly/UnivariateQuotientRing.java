package cc.redberry.rings.poly;

import cc.redberry.rings.AQuotientRing;
import cc.redberry.rings.WithVariables;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariateDivision;
import cc.redberry.rings.poly.univar.UnivariateDivision.InverseModMonomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;

/**
 * Univariate quotient ring
 */
public class UnivariateQuotientRing<Poly extends IUnivariatePolynomial<Poly>>
        extends AQuotientRing<Poly> implements IPolynomialRing<Poly> {
    private final Poly factory;
    public final Poly modulus;
    public final InverseModMonomial<Poly> fastDiv;

    public UnivariateQuotientRing(IPolynomialRing<Poly> baseRing, Poly modulus) {
        // fixme polys over Z ????
        super(baseRing);
        this.factory = modulus.createOne();
        this.modulus = modulus.clone();
        if (this.modulus.isOverField())
            this.modulus.monic();
        this.fastDiv = modulus.isOverField()
                ? UnivariateDivision.fastDivisionPreConditioning(this.modulus)
                : null;
    }

    @Override
    @SuppressWarnings("unchecked")
    public Poly mod(Poly el) {
        return (fastDiv == null
                ? modulus.isConstant()
                ? (Poly) divideByConstant((UnivariatePolynomial) el, (UnivariatePolynomial) modulus)
                : UnivariateDivision.pseudoDivideAndRemainder(el, modulus, true)[1]
                : UnivariateDivision.remainderFast(el, modulus, fastDiv, true))
                .canonical();
    }

    private static <E> UnivariatePolynomial<E> divideByConstant(UnivariatePolynomial<E> poly, UnivariatePolynomial<E> modulus) {
        assert modulus.isConstant();
        E cc = modulus.cc();
        return poly.mapCoefficients(poly.ring, cf -> poly.ring.remainder(cf, cc));
    }

    @Override
    public int nVariables() {
        return 1;
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
        if (variable != 0)
            throw new IllegalArgumentException();
        return factory.createMonomial(1);
    }

    @Override
    public String toString(String[] variables) {
        return toString(factory.coefficientRingToString(), variables);
    }

    @Override
    public String toString(String coefficientDomain, String[] variables) {
        return ((IPolynomialRing) baseRing).toString(coefficientDomain, variables) + "/<" + modulus.toString(variables) + ">";
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(nVariables()));
    }
}
