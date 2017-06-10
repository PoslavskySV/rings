package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCompositions;
import cc.r2.core.poly.CommonPolynomialsArithmetics;
import cc.r2.core.poly.univar.*;
import cc.r2.core.poly.univar.DivisionWithRemainder.*;
import cc.r2.core.util.ArraysUtil;

import java.util.*;

import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.asMultivariate;
import static cc.r2.core.poly.univar.DivisionWithRemainder.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class HenselLifting {
    private HenselLifting() {}


    private static <PolyZp extends IUnivariatePolynomial<PolyZp>>
    PolyZp[] solveDiophantine(PolyZp g, PolyZp f, PolyZp c) {
        PolyZp[] egcd = monicExtendedEuclid(g, f);
        PolyZp
                gcd = egcd[0],
                a = egcd[1],
                b = egcd[2],
                cGCD = divideExact(c, gcd, true),
                fGCD = divideExact(f, gcd, true),
                gGCD = divideExact(g, gcd, true);

        PolyZp[] qr = divideAndRemainder(a.clone().multiply(cGCD), fGCD, true);
        PolyZp q = qr[0], gSolution = qr[1];
        PolyZp fSolution = b.clone().multiply(cGCD).add(q.multiply(gGCD));

        assert g.clone().multiply(gSolution).add(f.clone().multiply(fSolution)).equals(c);
        return g.arrayNewInstance(gSolution, fSolution);
    }

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

    /* ======================================== Two-variable dense lifting ========================================== */


    /** Lift equation base = a * b mod x2 to actual solution */
    static void liftFactors(lMultivariatePolynomialZp base,
                            lMultivariatePolynomialZp a,
                            lMultivariatePolynomialZp b) {

        // a and b are coprime univariate polynomials over x1
        // we lift them up to the solution in (x1, x2)

        assert a.univariateVariable() == 0 && b.univariateVariable() == 0;
        assert base.nUsedVariables() == 2 && base.degree(0) > 0 && base.degree(1) > 0;
        assert base.evaluate(1, 0).equals(a.clone().multiply(b));

        lUnivariatePolynomialZp
                ua = a.asUnivariate(),
                ub = b.asUnivariate();
        UnivariatePolynomial<lMultivariatePolynomialZp>
                uBase = base.asUnivariate(1);

        assert AMultivariatePolynomial.asMultivariate(uBase, 1).equals(base);

        // solution of ua * s + ub * t = 1
        lUnivariatePolynomialZp
                egcd[] = monicExtendedEuclid(ua, ub),
                uaCofactor = egcd[1],
                ubCofactor = egcd[2];

        assert ua.clone().multiply(uaCofactor).add(ub.clone().multiply(ubCofactor)).isOne() : ua.clone().multiply(uaCofactor).add(ub.clone().multiply(ubCofactor));

        List<lUnivariatePolynomialZp>
                aCoefficients = new ArrayList<>(),
                bCoefficients = new ArrayList<>();
        aCoefficients.add(ua);
        bCoefficients.add(ub);

        InverseModMonomial<lUnivariatePolynomialZp>
                uaInvMod = fastDivisionPreConditioningWithLCCorrection(ua),
                ubInvMod = fastDivisionPreConditioningWithLCCorrection(ub);

        for (int degree = 1; degree <= base.degree(1); ++degree) {
            lUnivariatePolynomialZp rhs = uBase.get(degree).asUnivariate();
            for (int j = 1; j < degree; ++j)
                rhs = rhs.subtract(aCoefficients.get(j).clone().multiply(bCoefficients.get(degree - j)));

            if (rhs.isZero() && base.clone().subtract(a.clone().multiply(b)).isZero())
                break;

            lUnivariatePolynomialZp
                    aUpdate = ubCofactor.clone().multiply(rhs),
                    bUpdate = uaCofactor.clone().multiply(rhs);

            lUnivariatePolynomialZp[] qd = divideAndRemainderFast(aUpdate, ua, uaInvMod, false);
            aUpdate = qd[1];
            bUpdate = bUpdate.add(ub.clone().multiply(qd[0]));
//            aUpdate = remainderFast(aUpdate, ua, uaInvMod, false);
//            bUpdate = remainderFast(bUpdate, ub, ubInvMod, false);

            assert aUpdate.clone().multiply(ub).add(bUpdate.clone().multiply(ua)).equals(rhs);

            aCoefficients.add(aUpdate);
            bCoefficients.add(bUpdate);

            a.add(asMultivariate(aUpdate, 2, 0, a.ordering).multiply(a.createUnivariateMonomial(1, degree)));
            b.add(asMultivariate(bUpdate, 2, 0, b.ordering).multiply(b.createUnivariateMonomial(1, degree)));

            System.out.println("a  " + a.asUnivariate(0).lc());
            System.out.println("b  " + b.asUnivariate(0).lc());

            System.out.println(MultivariateReduction.dividesQ(b, b.asUnivariate(0).lc()));
        }
    }

    /** Lift equation base = a * b mod x2 to actual solution */
    static void liftFactors0(lMultivariatePolynomialZp base,
                             lMultivariatePolynomialZp a,
                             lMultivariatePolynomialZp b,
                             lMultivariatePolynomialZp aLC,
                             lMultivariatePolynomialZp bLC) {
        // a and b are coprime univariate polynomials over x1
        // we lift them up to the solution in (x1, x2)

        assert a.univariateVariable() == 0 && b.univariateVariable() == 0;
        assert base.nUsedVariables() == 2 && base.degree(0) > 0 && base.degree(1) > 0;
        assert base.evaluate(1, 0).equals(a.clone().multiply(b));

        lUnivariatePolynomialZp
                ua = a.asUnivariate(),
                ub = b.asUnivariate();

        if (!aLC.isOne() || !bLC.isOne()) {
            a.set(AMultivariatePolynomial.asMultivariate(a.asUnivariate(0).setLC(aLC), 0));
            b.set(AMultivariatePolynomial.asMultivariate(b.asUnivariate(0).setLC(bLC), 0));

            assert a.evaluate(1, 0).asUnivariate().equals(ua);
            assert b.evaluate(1, 0).asUnivariate().equals(ub);
            assert base.evaluate(1, 0).equals(a.clone().multiply(b).evaluate(1, 0));
        }

        // solution of ua * s + ub * t = 1
        lUnivariatePolynomialZp
                eGCD[] = monicExtendedEuclid(ua, ub),
                uaCoFactor = eGCD[1],
                ubCoFactor = eGCD[2];

        assert ua.clone().multiply(uaCoFactor).add(ub.clone().multiply(ubCoFactor)).isOne() : ua.clone().multiply(uaCoFactor).add(ub.clone().multiply(ubCoFactor));

        InverseModMonomial<lUnivariatePolynomialZp>
                uaInvMod = fastDivisionPreConditioningWithLCCorrection(ua),
                ubInvMod = fastDivisionPreConditioningWithLCCorrection(ub);

        for (int degree = 1; degree <= base.degree(1); ++degree) {
            lMultivariatePolynomialZp rhs = base.clone().subtract(a.clone().multiply(b));
            if (rhs.isZero())
                break;

            mod(rhs, 1, degree + 1);
            assert rhs.isZero() || (rhs.degrees(1).length == 1 && rhs.degrees(1)[0] == degree) : degree + " : " + rhs;

            lUnivariatePolynomialZp urhs = rhs.evaluate(1, 1).asUnivariate();

            lUnivariatePolynomialZp
                    aUpdate = ubCoFactor.clone().multiply(urhs),
                    bUpdate = uaCoFactor.clone().multiply(urhs);

            lUnivariatePolynomialZp[] qd = divideAndRemainderFast(aUpdate, ua, uaInvMod, false);
            aUpdate = qd[1];
            bUpdate = bUpdate.add(ub.clone().multiply(qd[0]));

            a.add(asMultivariate(aUpdate, 2, 0, a.ordering).multiply(a.createUnivariateMonomial(1, degree)));
            b.add(asMultivariate(bUpdate, 2, 0, b.ordering).multiply(b.createUnivariateMonomial(1, degree)));
        }
    }

    static void mod(lMultivariatePolynomialZp poly, int variable, int degree) {
        Iterator<Map.Entry<DegreeVector, lMonomialTerm>> it = poly.terms.entrySet().iterator();
        while (it.hasNext()) {
            lMonomialTerm term = it.next().getValue();
            if (term.exponents[variable] >= degree) {
                it.remove();
                poly.release();
            }
        }
    }

    /* =========================================== Single  lifting ============================================== */

    static final class Homomorphism {
        final int[] variables;
        final long[] values;

        Homomorphism(int[] variables, long[] values) {
            this.variables = variables;
            this.values = values;
        }

        long value(long variable) {
            for (int i = 0; i < variables.length; i++) {
                if (variables[i] == variable)
                    return values[i];
            }
            throw new RuntimeException("" + variable);
        }
    }

    static lUnivariatePolynomialZp cfx(
            lMultivariatePolynomialZp poly,
            Homomorphism morphism,
            int[] orders) {
        // todo : if morphism is zero => just pickup coefficients of with degree vectors = orders


//        if (orders[0] == 0 && orders[1] == 1 && orders[2] == 1 && orders[3] == 1)
//            System.out.println(poly);
        for (int i = 0; i < morphism.variables.length; ++i) {
            lMultivariatePolynomialZp tmp = poly;
            for (int j = 0; j < orders[i]; j++) {
                tmp = tmp.derivative(morphism.variables[i]);
            }
            poly = poly.derivative(morphism.variables[i], orders[i]);
            assert poly.equals(tmp);
        }
//        if (orders[0] == 0 && orders[1] == 1 && orders[2] == 1 && orders[3] == 1)
//            System.out.println(poly);

        lUnivariatePolynomialZp uPoly = poly.evaluate(morphism.variables, morphism.values).asUnivariate();
        for (int deg : orders)
            uPoly = uPoly.divide(uPoly.domain.factorial(deg));
        return uPoly;
    }


    public static lMultivariatePolynomialZp lift(
            UnivariatePolynomial<lMultivariatePolynomialZp> equation,
            Homomorphism morphism,
            lMultivariatePolynomialZp uSolution) {

        int nVariables = uSolution.nVariables;
        lMultivariatePolynomialZp interpolation = uSolution;
        lUnivariatePolynomialZp denom = equation.derivative().evaluate(uSolution).asUnivariate();
        for (int degree = 1; ; ++degree) {
            lMultivariatePolynomialZp eq = equation.evaluate(interpolation);
            if (eq.isZero())
                break;

            lMultivariatePolynomialZp update = uSolution.createZero();

            IntCompositions dPartitions = new IntCompositions(degree, morphism.variables.length);
            int[] dPartition;
//            System.out.println("======");
            while ((dPartition = dPartitions.take()) != null) {
//                System.out.println(Arrays.toString(dPartition));
                assert ArraysUtil.sum(dPartition) == degree;

                lUnivariatePolynomialZp num = cfx(eq, morphism, dPartition);
                lUnivariatePolynomialZp[] aCfx = DivisionWithRemainder.divideAndRemainder(num, denom, true);
                assert aCfx[1].isZero() : num + " / " + denom + " = { " + aCfx[0] + " ,  " + aCfx[1] + "}";

                lMultivariatePolynomialZp cfx = asMultivariate(aCfx[0].negate(), nVariables, 0, uSolution.ordering);
                out:
                for (int i = 0; i < morphism.variables.length; i++) {
                    int var = morphism.variables[i];
                    lMultivariatePolynomialZp diff = CommonPolynomialsArithmetics.polyPow(cfx.createLinear(var, cfx.domain.negate(morphism.value(var)), 1), dPartition[i], true);
                    cfx = cfx.multiply(diff);
                    if (cfx.isZero())
                        break;
                }
                update = update.add(cfx);
            }

//            System.out.println(interpolation);
            interpolation = interpolation.add(update);
        }

        return interpolation;
    }


    public static lMultivariatePolynomialZp[] liftFactorization(
            lMultivariatePolynomialZp base,
            Homomorphism morphism,
            lMultivariatePolynomialZp ua,
            lMultivariatePolynomialZp ub) {

        int nVariables = ua.nVariables;
        lMultivariatePolynomialZp
                aInterpolation = ua.clone(),
                bInterpolation = ub.clone();

        lUnivariatePolynomialZp
                ua0 = ua.asUnivariate(),
                ub0 = ub.asUnivariate();
        Comparator<DegreeVector> ordering = base.ordering;
        for (int degree = 1; ; ++degree) {
//            System.out.println("=== degree: " + degree);
            lMultivariatePolynomialZp diff = base.clone().subtract(aInterpolation.clone().multiply(bInterpolation));
            if (diff.isZero())
                break;
            lMultivariatePolynomialZp
                    aDelta = ua.createZero(),
                    bDelta = ub.createZero();

            IntCompositions dPartitions = new IntCompositions(degree, morphism.variables.length);
            int[] dPartition;
            while ((dPartition = dPartitions.take()) != null) {

                lUnivariatePolynomialZp rhs = cfx(diff, morphism, dPartition);
                lUnivariatePolynomialZp[] diophantine = solveDiophantine(ua0, ub0, rhs);
                lMultivariatePolynomialZp
                        aCfx = asMultivariate(diophantine[1], nVariables, 0, ordering),
                        bCfx = asMultivariate(diophantine[0], nVariables, 0, ordering);

                for (int i = 0; i < morphism.variables.length; i++) {
                    int var = morphism.variables[i];
                    lMultivariatePolynomialZp xDiff = CommonPolynomialsArithmetics.polyPow(base.createLinear(var, base.domain.negate(morphism.value(var)), 1), dPartition[i], true);
                    aCfx = aCfx.multiply(xDiff);
                    bCfx = bCfx.multiply(xDiff);
                    if (aCfx.isZero() && bCfx.isZero())
                        break;
                }

                aDelta.add(aCfx);
                bDelta.add(bCfx);
            }
            aInterpolation = aInterpolation.add(aDelta);
            bInterpolation = bInterpolation.add(bDelta);
        }

        return new lMultivariatePolynomialZp[]{aInterpolation, bInterpolation};
    }


}
