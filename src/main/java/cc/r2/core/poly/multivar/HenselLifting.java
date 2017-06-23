package cc.r2.core.poly.multivar;

import cc.r2.core.poly.Domain;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.MultivariatePolynomial.PrecomputedPowersHolder;
import cc.r2.core.poly.multivar.MultivariatePolynomial.USubstitution;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp.lPrecomputedPowersHolder;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp.lUSubstitution;
import cc.r2.core.poly.univar.*;
import cc.r2.core.poly.univar.DivisionWithRemainder.*;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;

import static cc.r2.core.poly.univar.DivisionWithRemainder.*;

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

    /** solves a * x + b * y = rhs */
    static <PolyZp extends IUnivariatePolynomial<PolyZp>>
    PolyZp[] solveDiophantine(PolyZp a, PolyZp b, PolyZp rhs) {
        PolyZp[] xgcd = monicExtendedEuclid(a, b);
        PolyZp
                x = xgcd[1].multiply(rhs),
                y = xgcd[2].multiply(rhs);

        PolyZp[] qd = divideAndRemainder(x, b, false);
        x = qd[1];
        y = y.add(qd[0].multiply(a));

//        PolyZp[] qd = divideAndRemainder(x, b, false);
//        x = qd[1];
//        qd = divideAndRemainder(y, a, false);
//        y = qd[1];

        qd[0] = x;
        qd[1] = y;
        return qd;
    }

    /* ======================================= 2-factor simple EZ lifting ========================================= */


    static void liftPair(lMultivariatePolynomialZp base,
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

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly primitivePart(Poly poly) {
        // multivariate GCDs will be used for calculation of primitive part
        return AMultivariatePolynomial.asMultivariate(poly.asUnivariate(0).primitivePart(), 0);
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

    /* ================================ 2-factor variable-by-variable EEZ lifting ================================== */

    /**
     * Holds a substitution x2 -> b2, ..., xN -> bN
     */
    interface IEvaluation<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {

        /**
         * Substitutes all variables starting from specified {@code variable} (inclusive), i.e.
         * {@code variable -> b_i, variable + 1-> b_(i + 1), ... xN -> bN}
         */
        Poly evaluateFrom(Poly poly, int variable);

        /**
         * @return {@code (1/order!) d(poly)/(d var) | var -> b}
         */
        Poly taylorCoefficient(Poly poly, int var, int order);

        /** @return (x_i - b_i)^exponent */
        Poly linearPower(int variable, int exponent);
    }

    static final class lEvaluation implements IEvaluation<lMonomialTerm, lMultivariatePolynomialZp> {
        final long[] values;
        final int nVariables;
        final lPrecomputedPowersHolder precomputedPowers;
        final lUSubstitution[] linearPowers;

        lEvaluation(int nVariables, long[] values, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.precomputedPowers = new lPrecomputedPowersHolder(nVariables, ArraysUtil.sequence(1, nVariables), values, domain);
            this.linearPowers = new lUSubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new lUSubstitution(lUnivariatePolynomialZ.create(-values[i], 1).modulus(domain), i + 1, nVariables, ordering);
        }

        @Override
        public lMultivariatePolynomialZp evaluateFrom(lMultivariatePolynomialZp poly, int variable) {
            return poly.evaluate(precomputedPowers, ArraysUtil.sequence(variable, nVariables));
        }

        @Override
        public lMultivariatePolynomialZp taylorCoefficient(lMultivariatePolynomialZp poly, int var, int order) {
            lMultivariatePolynomialZp derivative = poly.derivative(var, order);
            if (derivative.isZero())
                return derivative;
            return derivative
                    .evaluate(var, precomputedPowers.powers[var])
                    .divide(poly.domain.factorial(order));
        }

        @Override
        public lMultivariatePolynomialZp linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }
    }

    static final class Evaluation<E> implements IEvaluation<MonomialTerm<E>, MultivariatePolynomial<E>> {
        final E[] values;
        final int nVariables;
        final PrecomputedPowersHolder<E> precomputedPowers;
        final USubstitution<E>[] linearPowers;

        @SuppressWarnings("unchecked")
        Evaluation(int nVariables, E[] values, Domain<E> domain, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.precomputedPowers = new PrecomputedPowersHolder<>(nVariables, ArraysUtil.sequence(1, nVariables), values, domain);
            this.linearPowers = new USubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new USubstitution<>(UnivariatePolynomial.create(domain, domain.negate(values[i]), domain.getOne()), i + 1, nVariables, ordering);
        }

        @Override
        public MultivariatePolynomial<E> evaluateFrom(MultivariatePolynomial<E> poly, int variable) {
            return poly.evaluate(precomputedPowers, ArraysUtil.sequence(variable, nVariables));
        }

        @Override
        public MultivariatePolynomial<E> taylorCoefficient(MultivariatePolynomial<E> poly, int var, int order) {
            MultivariatePolynomial<E> derivative = poly.derivative(var, order);
            if (derivative.isZero())
                return derivative;
            return derivative
                    .evaluate(var, precomputedPowers.powers[var])
                    .divideOrNull(poly.domain.factorial(order));
        }

        @Override
        public MultivariatePolynomial<E> linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }
    }

    /** solves a * x + b * y = rhs for given univariate a, b and r (a and b are coprime) and unknown x and y */
    static final class UDiophantineSolver<uPoly extends IUnivariatePolynomial<uPoly>> {
        /** the given factors */
        final uPoly a, b;
        /** Bezout's factors: a * aCoFactor + b * bCoFactor = 1 */
        final uPoly aCoFactor, bCoFactor;

        UDiophantineSolver(uPoly a, uPoly b) {
            this.a = a;
            this.b = b;
            uPoly[] xgcd = monicExtendedEuclid(a, b);
            this.aCoFactor = xgcd[1];
            this.bCoFactor = xgcd[2];
        }

        /** the solution */
        uPoly x, y;

        void solve(uPoly rhs) {
            x = aCoFactor.clone().multiply(rhs);
            y = bCoFactor.clone().multiply(rhs);

            uPoly[] qd = DivisionWithRemainder.divideAndRemainder(x, b, false);
            x = qd[1];
            y = y.add(qd[0].multiply(a));
        }
    }

    /** solves a * x + b * y = rhs for given multivariate a, b and r (a and b are coprime) and unknown x and y */
    static final class DiophantineSolver<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        final IEvaluation<Term, Poly> evaluation;
        final Poly a, b;
        final Poly[] aImages, bImages;
        final UDiophantineSolver<uPoly> uSolver;
        final int[] degreeBounds;

        @SuppressWarnings("unchecked")
        public DiophantineSolver(Poly a,
                                 Poly b,
                                 IEvaluation<Term, Poly> evaluation,
                                 int[] degreeBounds) {
            this.a = a;
            this.b = b;
            this.evaluation = evaluation;
            this.degreeBounds = degreeBounds;

            aImages = a.arrayNewInstance(a.nVariables);
            bImages = a.arrayNewInstance(a.nVariables);

            for (int i = 0; i < a.nVariables; i++) {
                aImages[i] = evaluation.evaluateFrom(a, i + 1);
                bImages[i] = evaluation.evaluateFrom(b, i + 1);
            }

            uSolver = new UDiophantineSolver<>((uPoly) aImages[0].asUnivariate(), (uPoly) bImages[0].asUnivariate());
        }

        /** the solution */
        Poly x, y;

        @SuppressWarnings("unchecked")
        void solve(Poly rhs, int liftingVariable) {
            rhs = evaluation.evaluateFrom(rhs, liftingVariable + 1);
            if (liftingVariable == 0) {
                uSolver.solve((uPoly) rhs.asUnivariate());
                x = MultivariateGCD.asMultivariate(uSolver.x, rhs.nVariables, 0, rhs.ordering);
                y = MultivariateGCD.asMultivariate(uSolver.y, rhs.nVariables, 0, rhs.ordering);
                return;
            }

            // solve equation with x_i replaced with b_i:
            // a[x1, ..., x(i-1), b(i), ... b(N)] * x[x1, ..., x(i-1), b(i), ... b(N)]
            //    + b[x1, ..., x(i-1), b(i), ... b(N)] * y[x1, ..., x(i-1), b(i), ... b(N)]
            //         = rhs[x1, ..., x(i-1), b(i), ... b(N)]
            solve(rhs, liftingVariable - 1);

            // <- x and y are now:
            // x = x[x1, ..., x(i-1), b(i), ... b(N)]
            // y = y[x1, ..., x(i-1), b(i), ... b(N)]

            Poly
                    xSolution = x.clone(),
                    ySolution = y.clone();

            for (int degree = 1; degree <= degreeBounds[liftingVariable]; degree++) {
                // Δ = (rhs - a * x - b * y) mod (x_i - b_i)^degree
                Poly rhsDelta = rhs.clone().subtract(
                        aImages[liftingVariable].clone().multiply(xSolution)
                                .add(bImages[liftingVariable].clone().multiply(ySolution)));

                if (rhsDelta.isZero())
                    // we are done
                    break;

                rhsDelta = evaluation.taylorCoefficient(rhsDelta, liftingVariable, degree);

                solve(rhsDelta, liftingVariable - 1);
                assert x.isZero() || (x.degree(0) < b.degree(0)) : "\na:" + a + "\nb:" + b + "\nx:" + x + "\ny:" + y;

                // (x_i - b_i) ^ degree
                Poly idPower = evaluation.linearPower(liftingVariable, degree);
                xSolution.add(x.multiply(idPower));
                ySolution.add(y.multiply(idPower));
            }

            x = xSolution;
            y = ySolution;

            //assert assertSolution(rhs, liftingVariable);
        }

        boolean assertSolution(Poly rhs, int liftingVariable) {
            Poly delta = aImages[liftingVariable].clone().multiply(x).add(bImages[liftingVariable].clone().multiply(y)).subtract(rhs);
            for (int i = 0; i <= degreeBounds[liftingVariable]; i++)
                if (!evaluation.taylorCoefficient(delta, liftingVariable, i).isZero())
                    return false;
            return true;
        }
    }

    /**
     * Lifts factorization base = a * b mod <(x2-b2), ... , (xN-bN)>
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void liftWang(Poly base, Poly a, Poly b, IEvaluation<Term, Poly> evaluation) {
        Poly lc = base.lc(0);

        if (lc.isConstant())
            liftWang(base, a, b, null, null, evaluation);
        else {
            // imposing leading coefficients
            Poly lcCorrection = evaluation.evaluateFrom(lc, 1);

            assert !lcCorrection.isZero();
            assert a.lt().exponents[0] == a.degree(0);
            assert b.lt().exponents[0] == b.degree(0);

            a.monicWithLC(lcCorrection); // <- monic in x^n (depends on ordering!)
            b.monicWithLC(lcCorrection); // <- monic in x^n (depends on ordering!)
            base = base.clone().multiply(lc);

            liftWang(base, a, b, lc, lc, evaluation);

            a.set(primitivePart(a));
            b.set(primitivePart(b));
        }
    }

    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void liftWang(Poly base, Poly a, Poly b, Poly aLC, Poly bLC, IEvaluation<Term, Poly> evaluation) {
        // a and b are coprime univariate polynomials over x1
        // we lift them up to the solution in (x1, x2, ..., xN)
        assert a.univariateVariable() == 0 && b.univariateVariable() == 0 : a.univariateVariable() + "  " + b.univariateVariable();
        assert evaluation.evaluateFrom(base, 1).equals(evaluation.evaluateFrom(a.clone().multiply(b), 1))
                : "\n base: " + base + "\n a: " + a + "\n b: " + b;

        if (aLC != null) {
            // replace lc trick
            a.setLC(0, aLC);
        }

        if (bLC != null) {
            // replace lc trick
            b.setLC(0, bLC);
        }

        int[] degreeBounds = base.degrees();
        for (int liftingVariable = 1; liftingVariable < base.nVariables; ++liftingVariable) {
            DiophantineSolver<Term, Poly, ?> dSolver = new DiophantineSolver<>(
                    evaluation.evaluateFrom(a, liftingVariable),
                    evaluation.evaluateFrom(b, liftingVariable), evaluation, degreeBounds);

            // base[x1, x2, ..., x(i), b(i+1), ..., bN]
            Poly baseImage = evaluation.evaluateFrom(base, liftingVariable + 1);
            for (int degree = 1; degree <= degreeBounds[liftingVariable]; ++degree) {
                Poly rhsDelta =
                        baseImage.clone().subtract(a.clone().multiply(b.clone()));
                rhsDelta = evaluation.evaluateFrom(rhsDelta, liftingVariable + 1);
                if (rhsDelta.isZero())
                    break;
                rhsDelta = evaluation.taylorCoefficient(rhsDelta, liftingVariable, degree);
                assert rhsDelta.nUsedVariables() <= liftingVariable;

                dSolver.solve(rhsDelta, liftingVariable);

                // (x_i - b_i) ^ degree
                Poly idPower = evaluation.linearPower(liftingVariable, degree);
                a.add(dSolver.y.multiply(idPower));
                b.add(dSolver.x.multiply(idPower));
            }
        }
    }
}


