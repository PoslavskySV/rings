package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Ring;
import cc.redberry.rings.linear.LinearSolver;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64.lPrecomputedPowersHolder;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64.lUSubstitution;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.util.ArraysUtil;

import java.util.*;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.Frac;
import static cc.redberry.rings.Rings.PolynomialRing;

/**
 * Hensel lifting.
 *
 * @since 1.0
 */
public final class HenselLifting {
    private HenselLifting() {}

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static <PolyZp extends IUnivariatePolynomial<PolyZp>>
    PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
        PolyZp[] xgcd = UnivariateGCD.PolynomialExtendedGCD(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        assert xgcd[0].isConstant() : "bad xgcd: " + Arrays.toString(xgcd) + " for xgcd(" + a + ", " + b + ")";

        //normalize: x * a + y * b = 1
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[0].monic();

        return xgcd;
    }

    /**
     * Gives primitive part of poly considered as R[x2,x3,...,xN][x0]
     */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly primitivePart(Poly poly) {
        // multivariate GCDs will be used for calculation of primitive part
        return AMultivariatePolynomial.asMultivariate(poly.asUnivariate(0).primitivePart(), 0);
    }

    /**
     * Drops all terms of poly ∈ R[x1,x2,..,xN] which total degree with respect to [x2,.., xN] is equal or higher than
     * degree. NOTE: poly is not copied (returns the same reference)
     */
    static MultivariatePolynomialZp64 modImage(MultivariatePolynomialZp64 poly, int degree) {
        if (degree == 0)
            return poly.ccAsPoly();
        Iterator<Map.Entry<DegreeVector, MonomialZp64>> it = poly.terms.entrySet().iterator();
        while (it.hasNext()) {
            MonomialZp64 term = it.next().getValue();
            if (ArraysUtil.sum(term.exponents, 1) >= degree) {
                it.remove();
                poly.release();
            }
        }
        poly.release();
        return poly;
    }

    /**
     * Holds a substitution x2 -> b2, ..., xN -> bN
     */
    interface IEvaluation<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {

        /**
         * Substitutes all variables starting from specified {@code variable} (inclusive), i.e. {@code variable -> b_i,
         * variable + 1-> b_(i + 1), ... xN -> bN}
         */
        Poly evaluateFrom(Poly poly, int variable);

        /**
         * Substitutes all variables starting from specified variable {@code from} (inclusive) and except specified
         * variable
         */
        Poly evaluateFromExcept(Poly poly, int from, int except);

        /**
         * Substitute value for variable
         */
        Poly evaluate(Poly poly, int variable);

        /**
         * Sequentially evaliate all elements
         */
        default Poly[] evaluateFrom(Poly[] array, int variable) {
            Poly[] result = array[0].createArray(array.length);
            for (int i = 0; i < result.length; i++)
                result[i] = evaluateFrom(array[i], variable);
            return result;
        }

        boolean isZeroSubstitution(int variable);

        /**
         * @return {@code (1/order!) d(poly)/(d var) | var -> b}
         */
        default Poly taylorCoefficient(Poly poly, int variable, int order) {
            if (isZeroSubstitution(variable))
                return poly.coefficientOf(variable, order);

            return evaluate(poly.seriesCoefficient(variable, order), variable);
        }

        /** @return {@code (x_i - b_i)^exponent} */
        Poly linearPower(int variable, int exponent);

        default Poly modImage(Poly poly, int variable, int idealExponent) {
            if (idealExponent == 0)
                return poly.clone();
            int degree = poly.degree(variable);
            if (idealExponent < degree - idealExponent) {
                // select terms
                Poly result = poly.createZero();
                for (int i = 0; i < idealExponent; i++) {
                    Poly term = evaluate(poly.seriesCoefficient(variable, i), variable).multiply(linearPower(variable, i));
                    if (term.isZero())
                        continue;
                    result.add(term);
                }
                return result;
            } else {
                // drop terms
                poly = poly.clone();
                for (int i = idealExponent; i <= degree; i++) {
                    Poly term = evaluate(poly.seriesCoefficient(variable, i), variable).multiply(linearPower(variable, i));
                    if (term.isZero())
                        continue;
                    poly.subtract(term);
                }
                return poly;
            }
        }

        default Poly modImage(Poly poly, int[] degreeBounds) {
            for (int i = 1; i < degreeBounds.length; i++)
                poly = modImage(poly, i, degreeBounds[i] + 1);
            return poly;
        }

        IEvaluation<Term, Poly> dropVariable(int variable);

        IEvaluation<Term, Poly> renameVariables(int[] newVariablesExceptFirst);
    }

    /** Evaluations for {@link MultivariatePolynomialZp64} */
    static final class lEvaluation implements IEvaluation<MonomialZp64, MultivariatePolynomialZp64> {
        final long[] values;
        final int nVariables;
        final lPrecomputedPowersHolder precomputedPowers;
        final lUSubstitution[] linearPowers;
        final IntegersZp64 ring;
        final Comparator<DegreeVector> ordering;

        lEvaluation(int nVariables, long[] values, IntegersZp64 ring, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.ring = ring;
            this.ordering = ordering;
            this.precomputedPowers = MultivariatePolynomialZp64.mkPrecomputedPowers(nVariables, ring, ArraysUtil.sequence(1, nVariables), values);
            this.linearPowers = new lUSubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new lUSubstitution(UnivariatePolynomialZ64.create(-values[i], 1).modulus(ring), i + 1, nVariables, ordering);
        }

        @Override
        public MultivariatePolynomialZp64 evaluate(MultivariatePolynomialZp64 poly, int variable) {
            return poly.evaluate(variable, precomputedPowers.powers[variable]);
        }

        @Override
        public MultivariatePolynomialZp64 evaluateFrom(MultivariatePolynomialZp64 poly, int variable) {
            if (variable >= poly.nVariables)
                return poly.clone();
            if (variable == 1 && poly.univariateVariable() == 0)
                return poly.clone();
            return poly.evaluate(precomputedPowers, ArraysUtil.sequence(variable, nVariables));
        }

        @Override
        public MultivariatePolynomialZp64 evaluateFromExcept(MultivariatePolynomialZp64 poly, int from, int except) {
            if (from >= poly.nVariables)
                return poly.clone();
            if (from == 1 && poly.univariateVariable() == 0)
                return poly.clone();

            int[] vars = new int[poly.nVariables - from - 1];
            int c = 0;
            for (int i = from; i < except; i++)
                vars[c++] = i;
            for (int i = except + 1; i < nVariables; i++)
                vars[c++] = i;

            return poly.evaluate(precomputedPowers, vars);
        }

        @Override
        public MultivariatePolynomialZp64 linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }

        @Override
        public boolean isZeroSubstitution(int variable) {
            return values[variable - 1] == 0;
        }

        @Override
        public lEvaluation dropVariable(int variable) {
            return new lEvaluation(nVariables - 1, ArraysUtil.remove(values, variable - 1), ring, ordering);
        }

        @Override
        public lEvaluation renameVariables(int[] variablesExceptFirst) {
            return new lEvaluation(nVariables, map(values, variablesExceptFirst), ring, ordering);
        }

        @Override
        public String toString() {
            return Arrays.toString(values);
        }
    }

    private static long[] map(long[] oldArray, int[] mapping) {
        long[] newArray = new long[oldArray.length];
        for (int i = 0; i < oldArray.length; i++)
            newArray[i] = oldArray[mapping[i]];
        return newArray;
    }

    /** Generic evaluations */
    static final class Evaluation<E> implements IEvaluation<Monomial<E>, MultivariatePolynomial<E>> {
        final E[] values;
        final int nVariables;
        final Ring<E> ring;
        final MultivariatePolynomial.PrecomputedPowersHolder<E> precomputedPowers;
        final MultivariatePolynomial.USubstitution<E>[] linearPowers;
        final Comparator<DegreeVector> ordering;

        @SuppressWarnings("unchecked")
        Evaluation(int nVariables, E[] values, Ring<E> ring, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.ring = ring;
            this.ordering = ordering;
            this.precomputedPowers = MultivariatePolynomial.mkPrecomputedPowers(nVariables, ring, ArraysUtil.sequence(1, nVariables), values);
            this.linearPowers = new MultivariatePolynomial.USubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new MultivariatePolynomial.USubstitution<>(
                        UnivariatePolynomial.createUnsafe(ring, ring.createArray(ring.negate(values[i]), ring.getOne())),
                        i + 1, nVariables, ordering);
        }

        Evaluation<E> setRing(Ring<E> ring) {
            return new Evaluation<>(nVariables, values, ring, ordering);
        }

        @Override
        public MultivariatePolynomial<E> evaluate(MultivariatePolynomial<E> poly, int variable) {
            return poly.evaluate(variable, precomputedPowers.powers[variable]);
        }

        @Override
        public MultivariatePolynomial<E> evaluateFrom(MultivariatePolynomial<E> poly, int variable) {
            if (variable >= poly.nVariables)
                return poly.clone();
            if (variable == 1 && poly.univariateVariable() == 0)
                return poly.clone();
            return poly.evaluate(precomputedPowers, ArraysUtil.sequence(variable, nVariables));
        }

        @Override
        public MultivariatePolynomial<E> evaluateFromExcept(MultivariatePolynomial<E> poly, int from, int except) {
            if (from >= poly.nVariables)
                return poly.clone();
            if (from == 1 && poly.univariateVariable() == 0)
                return poly.clone();

            int[] vars = new int[poly.nVariables - from - 1];
            int c = 0;
            for (int i = from; i < except; i++)
                vars[c++] = i;
            for (int i = except + 1; i < nVariables; i++)
                vars[c++] = i;

            return poly.evaluate(precomputedPowers, vars);
        }

        @Override
        public MultivariatePolynomial<E> linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }

        @Override
        public boolean isZeroSubstitution(int variable) {
            return ring.isZero(values[variable - 1]);
        }

        @Override
        public Evaluation<E> dropVariable(int variable) {
            return new Evaluation<>(nVariables - 1, ArraysUtil.remove(values, variable - 1), ring, ordering);
        }

        @Override
        public Evaluation<E> renameVariables(int[] variablesExceptFirst) {
            return new Evaluation<>(nVariables, map(ring, values, variablesExceptFirst), ring, ordering);
        }

        @Override
        public String toString() {
            return Arrays.toString(values);
        }
    }

    private static <E> E[] map(Ring<E> ring, E[] oldArray, int[] mapping) {
        E[] newArray = ring.createArray(oldArray.length);
        for (int i = 0; i < oldArray.length; i++)
            newArray[i] = oldArray[mapping[i]];
        return newArray;
    }

    /**
     * Effectively holds all possible combinations of product of elements [p1,p2,...,pN]
     */
    static final class AllProductsCache<Poly extends IPolynomial<Poly>> {
        final Poly[] factors;
        final HashMap<BitSet, Poly> products = new HashMap<>();

        AllProductsCache(Poly[] factors) {
            assert factors.length >= 1;
            this.factors = factors;
        }

        private static BitSet clear(BitSet set, int from, int to) {
            set = (BitSet) set.clone();
            set.clear(from, to);
            return set;
        }

        Poly multiply(BitSet selector) {
            int cardinality = selector.cardinality();
            assert cardinality > 0;
            if (cardinality == 1)
                return factors[selector.nextSetBit(0)];
            Poly cached = products.get(selector);
            if (cached != null)
                return cached;
            // split BitSet into two ~equal parts:
            int half = cardinality / 2;
            for (int i = 0; ; ++i) {
                if (selector.get(i))
                    --half;
                if (half == 0) {
                    products.put(selector, cached =
                            multiply(clear(selector, 0, i + 1)).clone()
                                    .multiply(multiply(clear(selector, i + 1, factors.length))));
                    return cached;
                }
            }
        }

        Poly multiply(int[] selector) {
            BitSet bits = new BitSet(factors.length);
            for (int i = 0; i < selector.length; i++)
                bits.set(selector[i]);
            return multiply(bits);
        }

        int size() {
            return factors.length;
        }

        Poly get(int var) {
            return factors[var];
        }

        Poly except(int var) {
            BitSet bits = new BitSet(factors.length);
            bits.set(0, factors.length);
            bits.clear(var);
            return multiply(bits);
        }

        Poly from(int var) {
            BitSet bits = new BitSet(factors.length);
            bits.set(var, factors.length);
            return multiply(bits);
        }

        Poly[] exceptArray() {
            Poly[] arr = factors[0].createArray(factors.length);
            for (int i = 0; i < arr.length; i++)
                arr[i] = except(i);
            return arr;
        }

        Poly multiplyAll() {
            BitSet bits = new BitSet(factors.length);
            bits.set(0, factors.length);
            return multiply(bits);
        }
    }

    /**
     * Solves a * x + b * y = rhs for given univariate a, b and r (a and b are coprime) and unknown x and y
     */
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

            uPoly[] qd = UnivariateDivision.divideAndRemainder(x, b, false);
            x = qd[1];
            y = y.add(qd[0].multiply(a));
        }
    }

    /**
     * Solves a1 * x1 + a2 * x2 + ... = rhs for given univariate a1, a2, ... and rhs (all a_i are coprime) and unknown
     * x_i
     */
    static final class UMultiDiophantineSolver<uPoly extends IUnivariatePolynomial<uPoly>> {
        /** the given factors */
        final AllProductsCache<uPoly> factors;
        final UDiophantineSolver<uPoly>[] biSolvers;
        final uPoly[] solution;

        @SuppressWarnings("unchecked")
        UMultiDiophantineSolver(AllProductsCache<uPoly> factors) {
            this.factors = factors;
            this.biSolvers = new UDiophantineSolver[factors.size() - 1];
            for (int i = 0; i < biSolvers.length; i++)
                biSolvers[i] = new UDiophantineSolver<>(factors.get(i), factors.from(i + 1));
            this.solution = factors.factors[0].createArray(factors.factors.length);
        }

        void solve(uPoly rhs) {
            uPoly tmp = rhs.clone();
            for (int i = 0; i < factors.size() - 1; i++) {
                biSolvers[i].solve(tmp);
                solution[i] = biSolvers[i].y;
                tmp = biSolvers[i].x;
            }
            solution[factors.size() - 1] = tmp;
        }
    }

    /**
     * Solves a1 * x1 + a2 * x2 + ... = rhs for given multivariate a1, a2, ... and rhs (all a_i are coprime) and unknown
     * x_i
     */
    static final class MultiDiophantineSolver<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        final IEvaluation<Term, Poly> evaluation;
        final UMultiDiophantineSolver<uPoly> uSolver;
        final Poly[] solution;
        final int[] degreeBounds;
        final Ring<Poly> mRing;
        final UnivariatePolynomial<Poly>[][] imageSeries;

        MultiDiophantineSolver(IEvaluation<Term, Poly> evaluation,
                               Poly[] factors,
                               UMultiDiophantineSolver<uPoly> uSolver,
                               int[] degreeBounds,
                               int from) {
            assert from >= 1;
            this.evaluation = evaluation;
            this.uSolver = uSolver;
            this.degreeBounds = degreeBounds;
            Poly factory = factors[0];
            this.solution = factory.createArray(factors.length);
            this.mRing = new MultivariateRing<>(factors[0]);
            this.imageSeries = new UnivariatePolynomial[factory.nVariables][factors.length];
            for (int i = from - 1; i >= 1; --i)
                for (int j = 0; j < factors.length; j++)
                    this.imageSeries[i][j] = seriesExpansion(mRing, evaluation.evaluateFrom(factors[j], i + 1), i, evaluation);
        }

        void updateImageSeries(int liftingVariable, Poly[] factors) {
            for (int i = 0; i < factors.length; i++)
                imageSeries[liftingVariable][i] = seriesExpansion(mRing, factors[i], liftingVariable, evaluation);
        }

        @SuppressWarnings("unchecked")
        void solve(Poly rhs, int liftingVariable) {
            if (rhs.isZero()) {
                for (int i = 0; i < solution.length; i++)
                    solution[i] = rhs.createZero();
                return;
            }
            rhs = evaluation.evaluateFrom(rhs, liftingVariable + 1);
            if (liftingVariable == 0) {
                uSolver.solve((uPoly) rhs.asUnivariate());
                for (int i = 0; i < solution.length; i++)
                    solution[i] = AMultivariatePolynomial.asMultivariate(uSolver.solution[i], rhs.nVariables, 0, rhs.ordering);
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

            UnivariatePolynomial<Poly> rhsSeries = seriesExpansion(mRing, rhs, liftingVariable, evaluation);
            UnivariatePolynomial<Poly>[] tmpSolution = new UnivariatePolynomial[solution.length];
            for (int i = 0; i < tmpSolution.length; i++)
                tmpSolution[i] = seriesExpansion(mRing, solution[i], liftingVariable, evaluation);

            BernardinsTrick<Poly>[] pProducts = new BernardinsTrick[solution.length];
            for (int i = 0; i < solution.length; i++)
                pProducts[i] = createBernardinsTrick(new UnivariatePolynomial[]{tmpSolution[i], imageSeries[liftingVariable][i]}, degreeBounds[liftingVariable]);

            for (int degree = 1; degree <= degreeBounds[liftingVariable]; degree++) {
                // Δ = (rhs - a * x - b * y) mod (x_i - b_i)^degree
                Poly rhsDelta = rhsSeries.get(degree);
                for (int i = 0; i < solution.length; i++)
                    rhsDelta = rhsDelta.subtract(pProducts[i].fullProduct().get(degree));

                solve(rhsDelta, liftingVariable - 1);
                //assert x.isZero() || (x.degree(0) < b.degree(0)) : "\na:" + a + "\nb:" + b + "\nx:" + x + "\ny:" + y;

                // (x_i - b_i) ^ degree
                for (int i = 0; i < solution.length; i++)
                    pProducts[i].update(solution[i], rhs.createZero());
            }

            for (int i = 0; i < solution.length; i++)
                solution[i] = seriesToPoly(rhs, tmpSolution[i], liftingVariable, evaluation);
        }
    }

//    static <Poly extends IPolynomial<Poly>> void correctUnit(Poly poly, Poly[] factors) {
//        Poly lc = poly.lcAsPoly();
//        Poly flc = Arrays.stream(factors)
//                .map(IPolynomial::lcAsPoly)
//                .reduce(poly.createOne(), IPolynomial::multiply);
//        assert lc.isConstant();
//        assert flc.isConstant();
//
//        factors[0].multiplyByLC(lc.divideByLC(flc));
//    }

    /**
     * Fast bivariate Hensel lifting which uses dense representation for bivariate polynomials and without leading
     * coefficient correction (so the obtained factors may be not true factors but factorization mod I^bound holds)
     *
     * @param baseSeries  base bivariate polynomial in dense representation
     * @param factors     univariate factors
     * @param degreeBound bound on lifting degree
     * @return lifted bivariate factors in dense representation
     */
    @SuppressWarnings("unchecked")
    static <uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomial<uPoly>[] bivariateLiftDense(
            UnivariatePolynomial<uPoly> baseSeries, uPoly[] factors, int degreeBound) {
        AllProductsCache<uPoly> uFactors = new AllProductsCache<>(factors);
        // univariate multifactor diophantine solver
        UMultiDiophantineSolver<uPoly> uSolver = new UMultiDiophantineSolver<>(uFactors);

        UnivariatePolynomial<uPoly>[] solution = new UnivariatePolynomial[factors.length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = UnivariatePolynomial.constant(baseSeries.ring, uFactors.factors[i]);
            solution[i].ensureInternalCapacity(degreeBound + 1);
        }

        BernardinsTrickWithoutLCCorrection<uPoly> factorsProduct
                = new BernardinsTrickWithoutLCCorrection<>(solution);
        for (int degree = 1; degree <= degreeBound; ++degree) {
            uPoly rhsDelta = baseSeries.get(degree).clone().subtract(factorsProduct.fullProduct().get(degree));
            uSolver.solve(rhsDelta);
            factorsProduct.update(uSolver.solution);
        }

        return solution;
    }

    /**
     * Fast bivariate Hensel lifting which uses dense representation for bivariate polynomials
     *
     * @param base        the product
     * @param factors     univariate factors which will be lifted to true bivariate factors
     * @param evaluation  evaluation point
     * @param degreeBound bound on lifting degree
     */
    @SuppressWarnings("unchecked")
    public static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void bivariateLiftNoLCCorrection0(Poly base, Poly[] factors,
                                      IEvaluation<Term, Poly> evaluation, int degreeBound) {
        Ring<uPoly> uRing = new UnivariateRing<>((uPoly) factors[0].asUnivariate());
        UnivariatePolynomial<uPoly>[] res =
                bivariateLiftDense(seriesExpansionDense(uRing, base, 1, evaluation),
                        asUnivariate(factors, evaluation), degreeBound);
        for (int i = 0; i < res.length; i++)
            factors[i].set(denseSeriesToPoly(base, res[i], 1, evaluation));
    }

    /**
     * Impose leading coefficients on factors
     *
     * @param factors   the factors
     * @param factorsLC the correct leading coefficients
     */
    private static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void imposeLeadingCoefficients(Poly[] factors, Poly[] factorsLC) {
        if (factorsLC != null)
            for (int i = 0; i < factors.length; i++)
                if (factorsLC[i] != null)
                    factors[i].setLC(0, factorsLC[i]);
    }

    /**
     * Fast bivariate Hensel lifting which uses dense representation for bivariate polynomials
     *
     * @param base        the product
     * @param factors     univariate factors which will be lifted to true bivariate factors
     * @param factorsLC   leading coefficients (may be null)
     * @param evaluation  evaluation point
     * @param degreeBound bound on lifting degree
     */
    @SuppressWarnings("unchecked")
    static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void bivariateLift0(Poly base, Poly[] factors, Poly[] factorsLC,
                        IEvaluation<Term, Poly> evaluation, int degreeBound) {
        imposeLeadingCoefficients(factors, factorsLC);

        AllProductsCache<uPoly> uFactors = new AllProductsCache(asUnivariate(factors, evaluation));
        // univariate multifactor diophantine solver
        UMultiDiophantineSolver<uPoly> uSolver = new UMultiDiophantineSolver<>(uFactors);

        Ring<uPoly> uRing = new UnivariateRing<>(uFactors.get(0));
        UnivariatePolynomial<uPoly> baseSeries = seriesExpansionDense(uRing, base, 1, evaluation);
        UnivariatePolynomial<uPoly>[] solution = new UnivariatePolynomial[factors.length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = seriesExpansionDense(uRing, factors[i], 1, evaluation);
            solution[i].ensureInternalCapacity(degreeBound + 1);
        }

        BernardinsTrick<uPoly> product = createBernardinsTrick(solution, degreeBound);
        for (int degree = 1; degree <= degreeBound; ++degree) {
            uPoly rhsDelta = baseSeries.get(degree).clone().subtract(product.fullProduct().get(degree));
            uSolver.solve(rhsDelta);
            product.update(uSolver.solution);
        }

        for (int i = 0; i < solution.length; i++)
            factors[i].set(denseSeriesToPoly(base, solution[i], 1, evaluation));
    }


    @SuppressWarnings("unchecked")
    private static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    uPoly[] asUnivariate(Poly[] array, IEvaluation<Term, Poly> evaluation) {
        uPoly u0 = (uPoly) evaluation.evaluateFrom(array[0], 1).asUnivariate();
        uPoly[] res = u0.createArray(array.length);
        res[0] = u0;
        for (int i = 1; i < array.length; i++)
            res[i] = (uPoly) evaluation.evaluateFrom(array[i], 1).asUnivariate();
        return res;
    }

    @SuppressWarnings("unchecked")
    private static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    Poly[] asMultivariate(uPoly[] array, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        Poly u0 = AMultivariatePolynomial.asMultivariate(array[0], nVariables, variable, ordering);
        Poly[] res = u0.createArray(array.length);
        res[0] = u0;
        for (int i = 1; i < array.length; i++)
            res[i] = AMultivariatePolynomial.asMultivariate(array[i], nVariables, variable, ordering);
        return res;
    }

    /**
     * Generates a power series expansion for poly about the point specified by variable and evaluation
     */
    @SuppressWarnings("unchecked")
    public static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomial<uPoly> seriesExpansionDense(Ring<uPoly> ring, Poly poly, int variable, IEvaluation<Term, Poly> evaluate) {
        int degree = poly.degree(variable);
        uPoly[] coefficients = ring.createArray(degree + 1);
        for (int i = 0; i <= degree; i++)
            coefficients[i] = (uPoly) evaluate.taylorCoefficient(poly, variable, i).asUnivariate();
        return UnivariatePolynomial.createUnsafe(ring, coefficients);
    }

    /**
     * Converts power series {@link #seriesExpansionDense(Ring, AMultivariatePolynomial, int, IEvaluation)} back to
     * multivariate polynomial
     */
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    Poly denseSeriesToPoly(Poly factory, UnivariatePolynomial<uPoly> series, int seriesVariable, IEvaluation<Term, Poly> evaluation) {
        Poly result = factory.createZero();
        for (int i = 0; i <= series.degree(); i++) {
            Poly mPoly = AMultivariatePolynomial.asMultivariate(series.get(i), factory.nVariables, 0, factory.ordering);
            result = result.add(mPoly.multiply(evaluation.linearPower(seriesVariable, i)));
        }
        return result;
    }

    /**
     * Multivariate lift with automatic leading coefficient correction
     *
     * @param base       the product
     * @param factors    univariate factors which will be lifted to true bivariate factors
     * @param evaluation evaluation point
     */
    public static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void multivariateLiftAutomaticLC(Poly base,
                                     Poly[] factors,
                                     IEvaluation<Term, Poly> evaluation) {
        Poly lc = base.lc(0);

        if (lc.isConstant())
            multivariateLift0(base, factors, null, evaluation, base.degrees());
        else {
            // imposing leading coefficients
            Poly lcCorrection = evaluation.evaluateFrom(lc, 1);
            assert lcCorrection.isConstant();

            for (Poly factor : factors) {
                assert factor.lt().exponents[0] == factor.degree(0);
                factor.monicWithLC(lcCorrection.lcAsPoly());
            }

            base = base.clone().multiply(PolynomialMethods.polyPow(lc, factors.length - 1, true));

            multivariateLift0(base, factors, ArraysUtil.arrayOf(lc, factors.length), evaluation, base.degrees());

            for (Poly factor : factors)
                factor.set(primitivePart(factor));
        }
    }

    /**
     * Multivariate lift with automatic leading coefficient correction
     *
     * @param base       the product
     * @param factors    univariate factors which will be lifted to true bivariate factors
     * @param evaluation evaluation point
     */
    public static <
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void multivariateLiftAutomaticLC(Poly base,
                                     Poly[] factors,
                                     IEvaluation<Term, Poly> evaluation,
                                     int from) {
        Poly lc = base.lc(0);

        if (lc.isConstant())
            multivariateLift0(base, factors, null, evaluation, base.degrees());
        else {
            // imposing leading coefficients
            Poly lcCorrection = evaluation.evaluateFrom(lc, from);

            for (Poly factor : factors) {
                assert factor.lt().exponents[0] == factor.degree(0);
                factor.multiply(MultivariateDivision.divideExact(lcCorrection, factor.lc(0)));
            }

            base = base.clone().multiply(PolynomialMethods.polyPow(lc, factors.length - 1, true));

            multivariateLift0(base, factors, ArraysUtil.arrayOf(lc, factors.length), evaluation, base.degrees(), from);

            for (Poly factor : factors)
                factor.set(primitivePart(factor));
        }
    }


    /**
     * Multifactor multivariate Hensel lifting
     *
     * @param base         the product
     * @param factors      univariate factors which will be lifted to true bivariate factors
     * @param factorsLC    leading coefficients (may be null)
     * @param evaluation   evaluation point
     * @param degreeBounds bound on lifting degrees
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void multivariateLift0(Poly base,
                           Poly[] factors, Poly[] factorsLC,
                           IEvaluation<Term, Poly> evaluation,
                           int[] degreeBounds) {
        multivariateLift0(base, factors, factorsLC, evaluation, degreeBounds, 1);
    }

    /**
     * Sparsity threshold when to try sparse Hensel lifting
     */
    private static final double SPARSITY_THRESHOLD = 0.1;

    /**
     * Multifactor multivariate Hensel lifting
     *
     * @param base         the product
     * @param factors      univariate factors which will be lifted to true bivariate factors
     * @param factorsLC    leading coefficients (may be null)
     * @param evaluation   evaluation point
     * @param degreeBounds bound on lifting degrees
     * @param from         variable to start lifting with
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void multivariateLift0(Poly base,
                           Poly[] factors, Poly[] factorsLC,
                           IEvaluation<Term, Poly> evaluation,
                           int[] degreeBounds,
                           int from) {
        if (base.nVariables == 2) {
            bivariateLift0(base, factors, factorsLC, evaluation, degreeBounds[1]);
            return;
        }

        if (from == 2 && base.sparsity() < SPARSITY_THRESHOLD) {
            Poly[] sparseFactors = null;

            try {
                // ArithmeticException may raise in rare cases in Zp[X] where p is not a prime
                // (exception actually arises when calculating some GCDs, since some non-unit elements
                //  may be picked up from at random from Zp in randomized methods)
                sparseFactors = sparseLifting(base, factors, factorsLC);
            } catch (ArithmeticException ignore) { }

            if (sparseFactors != null) {
                System.arraycopy(sparseFactors, 0, factors, 0, factors.length);
                return;
            }
        }

        imposeLeadingCoefficients(factors, factorsLC);

        AllProductsCache<uPoly> uFactors = new AllProductsCache(asUnivariate(factors, evaluation));
        // univariate multifactor diophantine solver
        UMultiDiophantineSolver<uPoly> uSolver = new UMultiDiophantineSolver<>(uFactors);

        // initialize multivariate multifactor diophantine solver
        MultiDiophantineSolver<Term, Poly, ? extends IUnivariatePolynomial> dSolver = new MultiDiophantineSolver<>(
                evaluation,
                from == 1
                        ? (Poly[]) asMultivariate((IUnivariatePolynomial[]) uFactors.exceptArray(), base.nVariables, 0, base.ordering)
                        : new AllProductsCache<>(factors).exceptArray(),
                uSolver, degreeBounds, from);

        Ring<Poly> mRing = new MultivariateRing<>(factors[0]);
        for (int liftingVariable = from; liftingVariable < base.nVariables; ++liftingVariable) {
            Poly baseImage = evaluation.evaluateFrom(base, liftingVariable + 1);
            UnivariatePolynomial<Poly> baseSeries = seriesExpansion(mRing, baseImage, liftingVariable, evaluation);
            UnivariatePolynomial<Poly>[] solution = new UnivariatePolynomial[factors.length];
            for (int i = 0; i < solution.length; i++) {
                solution[i] = seriesExpansion(mRing, factors[i], liftingVariable, evaluation);
                solution[i].ensureInternalCapacity(degreeBounds[liftingVariable] + 1);
            }

            BernardinsTrick<Poly> product = createBernardinsTrick(solution, degreeBounds[liftingVariable]);
            for (int degree = 1; degree <= degreeBounds[liftingVariable]; ++degree) {
                Poly rhsDelta = baseSeries.get(degree).clone().subtract(product.fullProduct().get(degree));
                dSolver.solve(rhsDelta, liftingVariable - 1);
                product.update(dSolver.solution);
            }

            for (int i = 0; i < solution.length; i++)
                factors[i].set(seriesToPoly(base, solution[i], liftingVariable, evaluation));

            if (liftingVariable < base.nVariables) // don't perform on the last step
                dSolver.updateImageSeries(liftingVariable, new AllProductsCache<>(
                        evaluation.evaluateFrom(factors, liftingVariable + 1))
                        .exceptArray());
        }
    }

    /**
     * Generates a power series expansion for poly about the point specified by variable and evaluation
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    UnivariatePolynomial<Poly> seriesExpansion(Ring<Poly> ring, Poly poly, int variable, IEvaluation<Term, Poly> evaluate) {
        int degree = poly.degree(variable);
        Poly[] coefficients = ring.createArray(degree + 1);
        for (int i = 0; i <= degree; i++)
            coefficients[i] = evaluate.taylorCoefficient(poly, variable, i);
        return UnivariatePolynomial.createUnsafe(ring, coefficients);
    }

    /**
     * Converts power series {@link #seriesExpansion(Ring, AMultivariatePolynomial, int, IEvaluation)} back to
     * multivariate polynomial
     */
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly seriesToPoly(Poly factory, UnivariatePolynomial<Poly> series, int seriesVariable, IEvaluation<Term, Poly> evaluation) {
        Poly result = factory.createZero();
        for (int i = 0; i <= series.degree(); i++) {
            Poly mPoly = series.get(i);
            result = result.add(mPoly.multiply(evaluation.linearPower(seriesVariable, i)));
        }
        return result;
    }

    /*=========================== Bernardin's trick for fast error computation in lifting =============================*/

    static <Poly extends IPolynomial<Poly>>
    BernardinsTrick<Poly> createBernardinsTrick(UnivariatePolynomial<Poly>[] factors, int degreeBound) {
        if (Arrays.stream(factors).allMatch(UnivariatePolynomial::isConstant))
            return new BernardinsTrickWithoutLCCorrection<>(factors);
        else
            return new BernardinsTrickWithLCCorrection<>(factors, degreeBound);
    }

    static abstract class BernardinsTrick<Poly extends IPolynomial<Poly>> {
        // the factors
        final UnivariatePolynomial<Poly>[] factors;
        // partial products: {factor[0] * factor[1], (factor[0] * factor[1]) * factor[2], ... }
        final UnivariatePolynomial<Poly>[] partialProducts;
        final Ring<Poly> ring;

        BernardinsTrick(UnivariatePolynomial<Poly>... factors) {
            this.factors = factors;
            this.partialProducts = factors[0].createArray(factors.length - 1);
            this.ring = factors[0].ring;

            partialProducts[0] = factors[0].clone().multiply(factors[1]);
            for (int i = 1; i < partialProducts.length; i++)
                partialProducts[i] = partialProducts[i - 1].clone().multiply(factors[i + 1]);
        }

        /** update products */
        abstract void update(Poly... updates);

        /** get's (f_0 * f_1 * ... * f_N) mod y^(degree + 1) */
        final UnivariatePolynomial<Poly> fullProduct() {return partialProducts[partialProducts.length - 1];}

        final UnivariatePolynomial<Poly> partialProduct(int i) {
            return partialProducts[i];
        }
    }

    /** Bernardin's trick for fast f_0 * f_1 * ... * f_N computing (leading coefficients are discarded) */
    static final class BernardinsTrickWithoutLCCorrection<Poly extends IPolynomial<Poly>>
            extends BernardinsTrick<Poly> {
        public BernardinsTrickWithoutLCCorrection(UnivariatePolynomial<Poly>[] factors) {
            super(factors);
        }

        // current lift, so that factors are known mod I^pDegree
        private int pDegree = 0;

        /** update products */
        @Override
        void update(Poly... updates) {
            ++pDegree;
            // update factors
            for (int i = 0; i < factors.length; i++)
                factors[i].set(pDegree, updates[i]);

            // update the first product: factors[0] * factors[1]
            // k-th element is updated by (factors[0]_k * factors[1]_0 + factors[0]_0 * factors[1]_k)
            Poly updateValue = factors[0].get(pDegree).clone().multiply(factors[1].get(0))
                    .add(factors[1].get(pDegree).clone().multiply(factors[0].get(0)));
            partialProducts[0].addMonomial(updateValue, pDegree);

            // (k+1)-th element is calculated as (factors[0]_1*factors[1]_k + ... + factors[0]_k*factors[1]_1)
            Poly newElement = ring.getZero();
            for (int i = 1; i <= pDegree; i++)
                newElement.add(factors[0].get(i).clone().multiply(factors[1].get(pDegree - i + 1)));
            partialProducts[0].set(pDegree + 1, newElement);

            // => the first product (factors[0] * factors[1]) is updated
            // update other partial products accordingly
            for (int j = 1; j < partialProducts.length; j++) {

                // k-th element is updated by (update(p_k) * factors[j+1]_0 + p_0 * factors[j+1]_k),
                // where p is the previous partial product (without factors[j+1]) and
                // update(p_k) is the k-th element update of the previous partial product
                Poly currentUpdate =
                        partialProducts[j - 1].get(0).clone().multiply(factors[j + 1].get(pDegree))
                                .add(updateValue.multiply(factors[j + 1].get(0)));
                partialProducts[j].addMonomial(currentUpdate, pDegree);
                // cache current update for the next cycle
                updateValue = currentUpdate;

                // (k+1)-th element is calculated as (p[0]_1*factors[1]_k + ... + p[0]_k*factors[1]_1 + p[0]_(k+1)*factors[1]_0)
                newElement = ring.getZero();
                for (int i = 1; i <= (pDegree + 1); i++)
                    newElement.add(partialProducts[j - 1].get(i).clone().multiply(factors[j + 1].get(pDegree - i + 1)));
                partialProducts[j].set(pDegree + 1, newElement);
            }
        }
    }

    /** Bernardin's trick for fast f_0 * f_1 * ... * f_N computing (leading coefficients are took into account) */
    static final class BernardinsTrickWithLCCorrection<Poly extends IPolynomial<Poly>>
            extends BernardinsTrick<Poly> {
        final int degreeBound;

        BernardinsTrickWithLCCorrection(UnivariatePolynomial<Poly>[] factors, int degreeBound) {
            super(factors);
            this.degreeBound = degreeBound;
        }

        // current lift, so that factors are known mod I^pDegree
        private int pDegree = 0;

        private void updatePair(int iFactor, Poly leftUpdate, Poly rightUpdate, int degree) {
            // update factors
            UnivariatePolynomial<Poly>
                    left = iFactor == 0 ? factors[0] : partialProducts[iFactor - 1],
                    right = (iFactor + 1 < factors.length) ? factors[iFactor + 1] : null;

            if (leftUpdate != null) {
                left.addMonomial(leftUpdate, degree);
                if (iFactor < (factors.length - 1))
                    for (int i = degree; i <= Math.min(degreeBound, degree + right.degree()); i++)
                        updatePair(iFactor + 1, right.get(i - degree).clone().multiply(leftUpdate), null, i);
            }

            if (rightUpdate != null) {
                right.addMonomial(rightUpdate, degree);
                if (iFactor < (factors.length - 1))
                    for (int i = degree; i <= Math.min(degreeBound, degree + left.degree()); i++)
                        updatePair(iFactor + 1, left.get(i - degree).clone().multiply(rightUpdate), null, i);
            }
        }

        @Override
        void update(Poly... updates) {
            ++pDegree;
            updatePair(0, updates[0], updates[1], pDegree);
            for (int i = 0; i < factors.length - 2; i++)
                updatePair(i + 1, null, updates[i + 2], pDegree);
        }
    }


    /*=========================== Sparse Hensel lifting from bivariate factors =============================*/

    private static final long MAX_TERMS_IN_EXPAND_FORM = 8_388_608L;

    /**
     * Sparse Hensel lifting
     */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] sparseLifting(Poly base, Poly[] biFactors, Poly[] lc) {
        List<List<Term>>
                // terms that are fixed
                determinedTerms = new ArrayList<>(),
                // terms with unknowns
                undeterminedTerms = new ArrayList<>();

        // total number of unknown terms
        int nUnknowns = 0;
        long estimatedExpandedSize = 1;
        for (int i = 0; i < biFactors.length; ++i) {
            List<Term>
                    // fixed terms
                    fixed = new ArrayList<>(),
                    // unknown terms
                    unknown = new ArrayList<>();

            determinedTerms.add(fixed);
            undeterminedTerms.add(unknown);

            populateUnknownTerms(biFactors[i], lc == null ? base.createOne() : lc[i], fixed, unknown);
            nUnknowns += unknown.size();
            estimatedExpandedSize *= unknown.size();
        }

        if (estimatedExpandedSize > 1024L * base.size() || estimatedExpandedSize > MAX_TERMS_IN_EXPAND_FORM)
            // too large problem for sparse lifting -> probably we will run out of memory
            return null;

        // true factors represented as R[x0, x1, x2, ..., xN, u1, ..., uK]
        List<Poly> trueFactors = new ArrayList<>();
        int unkCounter = base.nVariables;
        for (int i = 0; i < biFactors.length; i++) {
            Poly trueFactor = base.createZero().joinNewVariables(nUnknowns);
            for (Term f : determinedTerms.get(i))
                trueFactor.add(f.joinNewVariables(nUnknowns));
            for (int j = 0; j < undeterminedTerms.get(i).size(); j++) {
                trueFactor.add(undeterminedTerms.get(i).get(j).joinNewVariables(nUnknowns)
                        .set(unkCounter, 1));
                ++unkCounter;
            }
            trueFactors.add(trueFactor);
        }

        // multiply our trueFactors in R[x0, x1, x2, ..., xN, u1, ..., uK]
        Poly lhsBase = trueFactors.stream().reduce(trueFactors.get(0).createOne(), (a, b) -> a.multiply(b));

        // <- matching lhsBase and base in (x0, x1)
        // base as R[x0, x1][x0, x1, x2, ... xN]
        MultivariatePolynomial<Poly> biBase =
                base.asOverMultivariate(ArraysUtil.sequence(2, base.nVariables));
        // The main equations to solve
        List<Equation<Term, Poly>> equations = new ArrayList<>();
        for (Monomial<Poly> rhs : biBase) {
            Poly cf = lhsBase.dropCoefficientOf(new int[]{0, 1}, Arrays.copyOf(rhs.exponents, 2));
            Equation<Term, Poly> eq = new Equation<>(cf, rhs.coefficient);

            if (!eq.isConsistent())
                // inconsistent equation -> bad lifting
                return null;

            if (!eq.isIdentity())
                // don't add identities
                equations.add(eq.canonical());
        }

        if (!lhsBase.isZero()) {
            Equation<Term, Poly> eq = new Equation<>(lhsBase, biBase.ring.getZero());
            if (!eq.isIdentity())
                equations.add(eq.canonical());
        }

        // all solutions we obtained so far
        ArrayList<BlockSolution<Term, Poly>> solutions = new ArrayList<>();
        main:
        while (!equations.isEmpty()) {
            // canonicalize all equations and rid out identical ones
            equations = new ArrayList<>(equations.stream().map(Equation::canonical).collect(Collectors.toSet()));

            // sort equations, so the first has less unknowns
            equations.sort(Comparator.comparingInt(eq -> eq.nUnknowns));

            assert equations.stream().allMatch(Equation::isConsistent);
            // filter linear equations
            List<Equation<Term, Poly>> linear =
                    equations.stream().filter(Equation::isLinear).collect(Collectors.toList());

            if (linear.isEmpty())
                // no any linear solutions more -> can't solve
                return null;

            // search for the block of linear equations
            List<Equation<Term, Poly>> block = new ArrayList<>();
            for (int i = 0; i < linear.size(); i++) {
                // take the base equation
                Equation<Term, Poly> baseEq = linear.get(i);
                block.add(baseEq);
                if (block.size() < baseEq.nUnknowns)
                    for (int j = 0; j < linear.size(); j++) {
                        if (i == j)
                            continue;

                        Equation<Term, Poly> eq2 = linear.get(j);
                        if (!eq2.hasOtherUnknownsThan(baseEq))
                            block.add(eq2);
                    }

                if (block.size() >= baseEq.nUnknowns) {
                    BlockSolution<Term, Poly> blockSolution = solveBlock(block);
                    if (blockSolution != null) {
                        solutions.add(blockSolution);
                        // remove all solved equations
                        equations.removeAll(block);
                        for (int k = 0; k < equations.size(); k++) {
                            Equation<Term, Poly> eq = equations.get(k).substituteSolutions(blockSolution);
                            if (!eq.isConsistent())
                                return null;
                            equations.set(k, eq);
                        }

                        assert equations.stream().allMatch(Equation::isConsistent);

                        // remove identity equations
                        equations.removeAll(equations.stream().filter(Equation::isIdentity).collect(Collectors.toList()));
                        continue main;
                    }
                }

                block.clear();
            }

            // no any solvable blocks
            return null;
        }

        return trueFactors.stream().map(
                factor -> {
                    // factor as R[x0, ..., xN][u1, ..., uK]
                    MultivariatePolynomial<Poly> result =
                            factor.asOverMultivariateEliminate(ArraysUtil.sequence(0, base.nVariables));
                    for (BlockSolution<Term, Poly> solution : solutions)
                        result = result.evaluate(solution.unknowns, solution.solutions);

                    assert result.isConstant();
                    return result.cc();
                }).toArray(base::createArray);
    }

    /**
     * Solves a block of linear equations in sparse lifting
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    BlockSolution<Term, Poly> solveBlock(List<Equation<Term, Poly>> block) {
        Equation<Term, Poly> baseEq = block.get(0);
        // unknown variables (indexed from zero)
        int[] unknowns = baseEq.getUnknowns();

        Poly factory = baseEq.rhs;
        // polynomial ring
        IPolynomialRing<Poly> polyRing = PolynomialRing(factory);
        // field of rational functions
        Ring<Rational<Poly>> fracRing = Frac(polyRing);

        // sort so the the first row will have maximal number of unknowns
        // this helps to solve system faster (less GCD's)
        block.sort(Comparator.comparing(eq -> -eq.nUnknowns));

        // lhs matrix
        List<Rational<Poly>[]> lhs = new ArrayList<>();
        // rhs column
        List<Rational<Poly>> rhs = new ArrayList<>();

        // try square system first
        for (int i = 0; i < unknowns.length; ++i)
            addRow(polyRing, block.get(i), lhs, rhs, unknowns);

        Rational<Poly>[] rSolution = new Rational[unknowns.length];
        while (true) {
            // convert to matrix
            Rational<Poly>[]
                    lhsMatrix[] = lhs.toArray(new Rational[lhs.size()][]),
                    rhsColumn = rhs.toArray(new Rational[rhs.size()]);

            LinearSolver.SystemInfo info = LinearSolver.solve(fracRing, lhsMatrix, rhsColumn, rSolution);
            if (info == LinearSolver.SystemInfo.Consistent)
                break;
            if (info == LinearSolver.SystemInfo.Inconsistent)
                return null;
            if (info == LinearSolver.SystemInfo.UnderDetermined) {
                // update matrix data with fresh reduced system
                lhs.clear();
                lhs.addAll(Arrays.asList(lhsMatrix));
                rhs.clear();
                rhs.addAll(Arrays.asList(rhsColumn));
                if (block.size() <= rhs.size())
                    return null;
                addRow(polyRing, block.get(rhs.size()), lhs, rhs, unknowns);
                continue;
            }
        }

        // real solution
        Poly[] solution = factory.createArray(rSolution.length);
        for (int i = 0; i < rSolution.length; i++) {
            Rational<Poly> r = rSolution[i];
            if (!r.denominator.isOne()) {
                // bad luck
                return null;
            }
            solution[i] = r.numerator;
        }

        BlockSolution<Term, Poly> blockSolution = new BlockSolution<>(unknowns, solution);
        // check the solution
        if (rhs.size() < block.size())
            for (int i = rhs.size(); i < block.size(); i++)
                if (!block.get(i).substituteSolutions(blockSolution).isConsistent())
                    return null;

        return blockSolution;
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void addRow(IPolynomialRing<Poly> polyRing, Equation<Term, Poly> eq,
                List<Rational<Poly>[]> lhs, List<Rational<Poly>> rhs, int[] unknowns) {
        Rational<Poly>[] row = new Rational[unknowns.length];
        for (int j = 0; j < row.length; j++) {
            // lhs matrix element
            MultivariatePolynomial<Poly> el = eq.lhs.coefficientOf(unknowns[j], 1);
            assert el.size() <= 1;
            row[j] = new Rational<>(polyRing, el.cc());
        }

        lhs.add(0, row);
        rhs.add(0, new Rational<>(polyRing, eq.rhs));
    }

    static final class BlockSolution<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        final int[] unknowns;
        final Poly[] solutions;

        BlockSolution(int[] unknowns, Poly[] solutions) {
            this.unknowns = unknowns;
            this.solutions = solutions;
        }

        @Override
        public String toString() {
            return "solution: " + Arrays.toString(solutions);

        }
    }

    /**
     * A single equation for sparse lifting
     */
    static final class Equation<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** lhs in R[x0, x1, ...., xN][u1, ... uL] */
        final MultivariatePolynomial<Poly> lhs;
        /** rhs in R[x0, x1, ...., xN] */
        final Poly rhs;
        /** lhs.degrees() */
        final int[] lhsDegrees;
        /** if equation is linear */
        final boolean isLinear;
        /** number of really unknown variables presented in this equation */
        final int nUnknowns;

        Equation(Poly lhs, Poly rhs) {
            this(lhs.asOverMultivariateEliminate(ArraysUtil.sequence(0, rhs.nVariables)), rhs);
        }

        Equation(MultivariatePolynomial<Poly> lhs, Poly rhs) {
            rhs.subtract(lhs.cc());
            lhs.subtract(lhs.cc());
            this.lhs = lhs;
            this.rhs = rhs;
            this.lhsDegrees = this.lhs.degrees();
            this.isLinear = lhs.getSkeleton().stream().allMatch(__ -> __.totalDegree <= 1);
            int nUnknowns = 0;
            for (int lhsDegree : lhsDegrees)
                nUnknowns += lhsDegree == 0 ? 0 : 1;
            this.nUnknowns = nUnknowns;
        }

        boolean isConsistent() {
            return !isIdentity() || lhs.cc().equals(rhs);
        }

        boolean hasOtherUnknownsThan(Equation<Term, Poly> other) {
            for (int i = 0; i < lhsDegrees.length; i++)
                if (lhsDegrees[i] > other.lhsDegrees[i])
                    return true;
            return false;
        }

        boolean isIdentity() { return ArraysUtil.sum(lhsDegrees) == 0;}

        boolean isLinear() { return isLinear;}

        int[] getUnknowns() {
            int[] unknowns = new int[nUnknowns];
            int i = -1;
            for (int j = 0; j < lhsDegrees.length; j++)
                if (lhsDegrees[j] != 0)
                    unknowns[++i] = j;
            return unknowns;
        }

        Equation<Term, Poly> substituteSolutions(BlockSolution<Term, Poly> solution) {
            // IMPORTANT: invocation of clone() is very important since rhs is modified in constructor!
            return new Equation<>(lhs.evaluate(solution.unknowns, solution.solutions), rhs.clone());
        }

        Equation<Term, Poly> canonical() {
            Poly gcd = MultivariateGCD.PolynomialGCD(lhs.content(), rhs);
            if (rhs.isOverField() && !rhs.isZero()) {
                lhs.divideExact(rhs.lcAsPoly());
                rhs.monic();
            }
            return new Equation<>(lhs.divideExact(gcd), MultivariateDivision.divideExact(rhs, gcd));
        }

        @Override
        public String toString() {
            return lhs.toString() + " = " + rhs;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Equation<?, ?> equation = (Equation<?, ?>) o;

            if (!lhs.equals(equation.lhs)) return false;
            return rhs.equals(equation.rhs);
        }

        @Override
        public int hashCode() {
            int result = lhs.hashCode();
            result = 31 * result + rhs.hashCode();
            return result;
        }
    }

    /** split terms in polynomials that are fixed (those coming from l.c.) and with unknown coefficients */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void populateUnknownTerms(Poly biPoly, Poly lc, List<Term> fixed, List<Term> unknown) {
        // degree in x0
        int xDeg = biPoly.degree(0);
        for (Term term : biPoly) {
            if (term.exponents[0] == xDeg) {
                Poly cf = lc.coefficientOf(1, term.exponents[1]);
                term = term.setCoefficientFrom(cf.monomialAlgebra.getUnitTerm(cf.nVariables));
                fixed.addAll(cf.multiply(term).collection());
            } else
                unknown.add(term);
        }
    }
}


