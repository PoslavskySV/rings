package cc.redberry.rings.scaladsl;

import cc.redberry.rings.poly.FiniteField;
import cc.redberry.rings.poly.multivar.AMonomial;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateDivision;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import scala.NotImplementedError;

/** Ugly type erasures for interoperability with Java */
final class JavaConversions {
    private JavaConversions() {}

    @SuppressWarnings("unchecked")
    static <Poly> Poly[] divideAndRemainder(Poly dividend, Poly... dividers) {
        return (Poly[]) MultivariateDivision.divideAndRemainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) dividers);
    }

    @SuppressWarnings("unchecked")
    static <Poly> Poly remainder(Poly dividend, Poly... dividers) {
        return (Poly) MultivariateDivision.remainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) dividers);
    }

    @SuppressWarnings("unchecked")
    static <Poly> Poly swapVariables(Poly poly, int i, int j) {
        return (Poly) AMultivariatePolynomial.swapVariables((AMultivariatePolynomial) poly, i, j);
    }

    @SuppressWarnings("unchecked")
    static <E extends IUnivariatePolynomial<E>, C> SimpleFieldExtension<E, C>
    mkScalaFieldExtension(cc.redberry.rings.poly.IPolynomialRing<E> javaRing, String variable) {
        if (javaRing instanceof FiniteField && javaRing.factory() instanceof UnivariatePolynomialZp64)
            return (SimpleFieldExtension<E, C>) new GaloisField64((FiniteField<UnivariatePolynomialZp64>) javaRing, variable);
        if (javaRing instanceof FiniteField)
            return (SimpleFieldExtension<E, C>) new GaloisField((FiniteField) javaRing, variable, null);
        if (javaRing instanceof cc.redberry.rings.poly.AlgebraicNumberField)
            return (SimpleFieldExtension<E, C>) new AlgebraicNumberField((cc.redberry.rings.poly.AlgebraicNumberField) javaRing, variable, null);

        throw new NotImplementedError();
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            mPoly extends AMultivariatePolynomial<Term, mPoly>,
            sPoly extends IUnivariatePolynomial<sPoly>,
            C> MultipleFieldExtension<Term, mPoly, sPoly, C>
    mkSimpleToMultiple(SimpleFieldExtension<sPoly, C> ext) {
        IPolynomialRing pr = (IPolynomialRing) ext;
        MultipleFieldExtension<Term, mPoly, sPoly, C> result = new MultipleFieldExtension((cc.redberry.rings.poly.MultipleFieldExtension)
                ext.implicitConversions().asMultipleExtension(), pr.variables());
        result.coder().withEncoder(pr.coder());
        return result;
    }
}
