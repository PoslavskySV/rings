package cc.r2.core.poly.multivar;

import cc.r2.core.poly.UnivariatePolynomials;
import cc.r2.core.poly.univar.DivisionWithRemainder.InverseModMonomial;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariateGCD;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;

import static cc.r2.core.poly.univar.DivisionWithRemainder.divideAndRemainderFast;
import static cc.r2.core.poly.univar.DivisionWithRemainder.fastDivisionPreConditioningWithLCCorrection;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class HenselLifting {
    private HenselLifting() {}


    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static <PolyZp extends IUnivariatePolynomial<PolyZp>>
    PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
        PolyZp[] xgcd = UnivariateGCD.ExtendedEuclid(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        assert xgcd[0].isConstant() : "bad xgcd: " + Arrays.toString(xgcd) + " for xgcd(" + a + ", " + b + ")";

        //normalize: x * a + y * b = 1
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[0].monic();

        return xgcd;
    }

    /* ======================================== 2-factor lifting ========================================== */


    public static void liftPair(lMultivariatePolynomialZp base,
                                lMultivariatePolynomialZp a,
                                lMultivariatePolynomialZp b) {
        lMultivariatePolynomialZp lc = base.lc(0);

        if (lc.isConstant())
            liftPair0(base, a, b, null, null);
        else {
            // imposing leading coefficients
            lMultivariatePolynomialZp lcCorrection = lc.ccAsPoly();

            assert a.lt().exponents[0] == a.degree(0);
            assert b.lt().exponents[0] == b.degree(0);

            a.monic(lcCorrection.lc()); // <- monic in x^n (depends on ordering!)
            b.monic(lcCorrection.lc()); // <- monic in x^n (depends on ordering!)
            base = base.clone().multiply(lc);

            liftPair0(base, a, b, lc, lc);

            a.set(primitivePart(a));
            b.set(primitivePart(b));
        }
    }

     static lMultivariatePolynomialZp primitivePart(lMultivariatePolynomialZp poly) {
        if (poly.nVariables == 2)
            // univariate GCDs will be used for calculation of primitive part
            return toSparseRepresentation(poly, toDenseRepresentation(poly).primitivePart());
        else
            // multivariate GCDs will be used for calculation of primitive part
            return AMultivariatePolynomial.asMultivariate(poly.asUnivariate(0).primitivePart(), 0);
    }

    /** dense representation R[x][y] of bivariate polynomial **/
    private static UnivariatePolynomial<lUnivariatePolynomialZp> toDenseRepresentation(lMultivariatePolynomialZp poly) {
        UnivariatePolynomials<lUnivariatePolynomialZp> domain = new UnivariatePolynomials<>(lUnivariatePolynomialZp.zero(poly.domain));
        lUnivariatePolynomialZp[] univarData = domain.createZeroesArray(poly.degree(0) + 1);
        for (lMonomialTerm e : poly.terms)
            univarData[e.exponents[0]].addMonomial(e.coefficient, e.exponents[1]);
        return UnivariatePolynomial.createUnsafe(domain, univarData);
    }

    /** sparse representation R[x, y] of bivariate polynomial R[x][y] **/
    private static lMultivariatePolynomialZp toSparseRepresentation(lMultivariatePolynomialZp factory, UnivariatePolynomial<lUnivariatePolynomialZp> poly) {
        lMultivariatePolynomialZp r = factory.createZero();
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;

            lUnivariatePolynomialZp yPoly = poly.get(i);
            r.add(lMultivariatePolynomialZp.asMultivariate(yPoly, 2, 1, factory.ordering).multiplyByMonomial(0, i));
        }
        return r;
    }

    /** Lift equation base = a * b mod x2 to actual solution */
    static void liftPair0(lMultivariatePolynomialZp base,
                          lMultivariatePolynomialZp a,
                          lMultivariatePolynomialZp b,
                          lMultivariatePolynomialZp aLC,
                          lMultivariatePolynomialZp bLC) {
        // a and b are coprime univariate polynomials over x1
        // we lift them up to the solution in (x1, x2)

        assert a.univariateVariable() == 0 && b.univariateVariable() == 0 : a.univariateVariable() + "  " + b.univariateVariable();
        assert modImage(base.clone(), 1).equals(a.clone().multiply(b));

        lUnivariatePolynomialZp
                ua = a.asUnivariate(),
                ub = b.asUnivariate();

        if (aLC != null) {
            // replace lc trick
            a.setLC(0, aLC);
            assert modImage(a.clone(), 1).asUnivariate().equals(ua);
        }

        if (bLC != null) {
            // replace lc trick
            b.setLC(0, bLC);
            assert modImage(b.clone(), 1).asUnivariate().equals(ub);
        }

        assert modImage(base.clone(), 1).equals(modImage(a.clone().multiply(b), 1));

        // solution of ua * s + ub * t = 1
        lUnivariatePolynomialZp
                eGCD[] = monicExtendedEuclid(ua, ub),
                uaCoFactor = eGCD[1],
                ubCoFactor = eGCD[2];

        InverseModMonomial<lUnivariatePolynomialZp>
                uaInvMod = fastDivisionPreConditioningWithLCCorrection(ua);

        int maxDegree = ArraysUtil.sum(base.degrees(), 1);
        for (int degree = 1; degree <= maxDegree; ++degree) {
            // reduce a and b mod degree to make things faster
            lMultivariatePolynomialZp
                    aMod = a.clone(),
                    bMod = b.clone();
            modImage(aMod, degree + 1);
            modImage(bMod, degree + 1);

            lMultivariatePolynomialZp rhs = base.clone().subtract(aMod.multiply(bMod));
            if (rhs.isZero())
                break;

            modImage(rhs, degree + 1);
            MultivariatePolynomial<lUnivariatePolynomialZp> rhsMod = rhs.asOverUnivariate(0);
            for (MonomialTerm<lUnivariatePolynomialZp> term : rhsMod) {
                lUnivariatePolynomialZp urhs = term.coefficient;

                lUnivariatePolynomialZp
                        aUpdate = ubCoFactor.clone().multiply(urhs),
                        bUpdate = uaCoFactor.clone().multiply(urhs);

                lUnivariatePolynomialZp[] qd = divideAndRemainderFast(aUpdate, ua, uaInvMod, false);
                aUpdate = qd[1];
                bUpdate = bUpdate.add(qd[0].multiply(ub));

                a.add(lMultivariatePolynomialZp
                        .asMultivariate(aUpdate, base.nVariables, 0, base.ordering)
                        .multiplyByDegreeVector(term));

                b.add(lMultivariatePolynomialZp
                        .asMultivariate(bUpdate, base.nVariables, 0, base.ordering)
                        .multiplyByDegreeVector(term));
            }
        }
    }

    /**
     * Drops all terms of poly ∈ R[x1,x2,..,xN] which total degree with respect to [x2,.., xN] is equal or higher
     * than degree. NOTE: poly is not copied (returns the same reference)
     */
    static lMultivariatePolynomialZp modImage(lMultivariatePolynomialZp poly, int degree) {
        if (degree == 0)
            return poly.ccAsPoly();
        Iterator<Map.Entry<DegreeVector, lMonomialTerm>> it = poly.terms.entrySet().iterator();
        while (it.hasNext()) {
            lMonomialTerm term = it.next().getValue();
            if (ArraysUtil.sum(term.exponents, 1) >= degree) {
                it.remove();
                poly.release();
            }
        }
        return poly;
    }
}
