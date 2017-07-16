package cc.r2.core.poly.multivar;

import cc.r2.core.poly.*;
import cc.r2.core.poly.multivar.MultivariatePolynomial.PrecomputedPowersHolder;
import cc.r2.core.poly.multivar.MultivariatePolynomial.USubstitution;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp.lPrecomputedPowersHolder;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp.lUSubstitution;
import cc.r2.core.poly.univar.*;
import cc.r2.core.util.ArraysUtil;

import java.util.*;

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

    /**
     * Gives primitive part of poly considered as R[x2,x3,...,xN][x0]
     */
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly primitivePart(Poly poly) {
        // multivariate GCDs will be used for calculation of primitive part
        return AMultivariatePolynomial.asMultivariate(poly.asUnivariate(0).primitivePart(), 0);
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
         * Substitute value for variable
         */
        Poly evaluate(Poly poly, int variable);

        /**
         * Sequentially evaliate all elements
         */
        default Poly[] evaluateFrom(Poly[] array, int variable) {
            Poly[] result = array[0].arrayNewInstance(array.length);
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

        /** @return (x_i - b_i)^exponent */
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
    }

    /** Evaluations for {@link lMultivariatePolynomialZp} */
    static final class lEvaluation implements IEvaluation<lMonomialTerm, lMultivariatePolynomialZp> {
        final long[] values;
        final int nVariables;
        final lPrecomputedPowersHolder precomputedPowers;
        final lUSubstitution[] linearPowers;
        final lIntegersModulo domain;
        final Comparator<DegreeVector> ordering;

        lEvaluation(int nVariables, long[] values, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.domain = domain;
            this.ordering = ordering;
            this.precomputedPowers = new lPrecomputedPowersHolder(nVariables, ArraysUtil.sequence(1, nVariables), values, domain);
            this.linearPowers = new lUSubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new lUSubstitution(lUnivariatePolynomialZ.create(-values[i], 1).modulus(domain), i + 1, nVariables, ordering);
        }

        @Override
        public lMultivariatePolynomialZp evaluate(lMultivariatePolynomialZp poly, int variable) {
            return poly.evaluate(variable, precomputedPowers.powers[variable]);
        }

        @Override
        public lMultivariatePolynomialZp evaluateFrom(lMultivariatePolynomialZp poly, int variable) {
            if (variable >= poly.nVariables)
                return poly.clone();
            if (variable == 1 && poly.univariateVariable() == 0)
                return poly.clone();
            return poly.evaluate(precomputedPowers, ArraysUtil.sequence(variable, nVariables));
        }

        @Override
        public lMultivariatePolynomialZp linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }

        @Override
        public boolean isZeroSubstitution(int variable) {
            return values[variable - 1] == 0;
        }

        @Override
        public lEvaluation dropVariable(int variable) {
            return new lEvaluation(nVariables - 1, ArraysUtil.remove(values, variable - 1), domain, ordering);
        }
    }

    /** Generic evaluations */
    static final class Evaluation<E> implements IEvaluation<MonomialTerm<E>, MultivariatePolynomial<E>> {
        final E[] values;
        final int nVariables;
        final Domain<E> domain;
        final PrecomputedPowersHolder<E> precomputedPowers;
        final USubstitution<E>[] linearPowers;
        final Comparator<DegreeVector> ordering;

        @SuppressWarnings("unchecked")
        Evaluation(int nVariables, E[] values, Domain<E> domain, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.values = values;
            this.domain = domain;
            this.ordering = ordering;
            this.precomputedPowers = new PrecomputedPowersHolder<>(nVariables, ArraysUtil.sequence(1, nVariables), values, domain);
            this.linearPowers = new USubstitution[nVariables - 1];
            for (int i = 0; i < nVariables - 1; i++)
                linearPowers[i] = new USubstitution<>(UnivariatePolynomial.create(domain, domain.negate(values[i]), domain.getOne()), i + 1, nVariables, ordering);
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
        public MultivariatePolynomial<E> linearPower(int variable, int exponent) {
            return linearPowers[variable - 1].pow(exponent);
        }

        @Override
        public boolean isZeroSubstitution(int variable) {
            return domain.isZero(values[variable - 1]);
        }

        @Override
        public Evaluation<E> dropVariable(int variable) {
            return new Evaluation<>(nVariables - 1, ArraysUtil.remove(values, variable - 1), domain, ordering);
        }
    }

    /**
     * Effectively holds all possible combinations of product of elements [p1,p2,...,pN]
     */
    static final class AllProductsCache<Poly extends IGeneralPolynomial<Poly>> {
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
            Poly[] arr = factors[0].arrayNewInstance(factors.length);
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

            uPoly[] qd = DivisionWithRemainder.divideAndRemainder(x, b, false);
            x = qd[1];
            y = y.add(qd[0].multiply(a));
        }
    }

    /**
     * Solves a1 * x1 + a2 * x2 + ... = rhs for given univariate a1, a2, ... and rhs (all a_i are coprime) and unknown x_i
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
            this.solution = factors.factors[0].arrayNewInstance(factors.factors.length);
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
     * Solves a1 * x1 + a2 * x2 + ... = rhs for given multivariate a1, a2, ... and rhs (all a_i are coprime) and unknown x_i
     */
    static final class MultiDiophantineSolver<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        final IEvaluation<Term, Poly> evaluation;
        final UMultiDiophantineSolver<uPoly> uSolver;
        final Poly[] solution;
        final Poly[][] images;
        final int[] degreeBounds;

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
            this.solution = factory.arrayNewInstance(factors.length);
            this.images = factory.arrayNewInstance2D(factory.nVariables, factors.length);
            this.images[from - 1] = factors;
            for (int i = from - 2; i >= 0; --i)
                this.images[i] = evaluation.evaluateFrom(this.images[i + 1], i + 1);
        }

        @SuppressWarnings("unchecked")
        void solve(Poly rhs, int liftingVariable) {
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

            Poly[] tmpSolution = solution[0].arrayNewInstance(solution.length);
            for (int i = 0; i < tmpSolution.length; i++)
                tmpSolution[i] = solution[i].clone();

            for (int degree = 1; degree <= degreeBounds[liftingVariable]; degree++) {
                // Δ = (rhs - a * x - b * y) mod (x_i - b_i)^degree
                Poly rhsDelta = rhs.clone();
                for (int i = 0; i < solution.length; i++)
                    rhsDelta = rhsDelta.subtract(images[liftingVariable][i].clone().multiply(tmpSolution[i]));

                if (rhsDelta.isZero())
                    // we are done
                    break;

                rhsDelta = evaluation.taylorCoefficient(rhsDelta, liftingVariable, degree);

                solve(rhsDelta, liftingVariable - 1);
                //assert x.isZero() || (x.degree(0) < b.degree(0)) : "\na:" + a + "\nb:" + b + "\nx:" + x + "\ny:" + y;

                // (x_i - b_i) ^ degree
                Poly idPower = evaluation.linearPower(liftingVariable, degree);
                for (int i = 0; i < tmpSolution.length; i++)
                    tmpSolution[i].add(solution[i].multiply(idPower));
            }

            System.arraycopy(tmpSolution, 0, solution, 0, tmpSolution.length);
        }
    }

    static <Poly extends IGeneralPolynomial<Poly>> void correctUnit(Poly poly, Poly[] factors) {
        Poly lc = poly.lcAsPoly();
        Poly flc = Arrays.stream(factors)
                .map(IGeneralPolynomial::lcAsPoly)
                .reduce(poly.createOne(), IGeneralPolynomial::multiply);
        assert lc.isConstant();
        assert flc.isConstant();

        factors[0].multiplyByLC(lc.divideByLC(flc));
    }

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
        AllProductsCache<uPoly> uFactors = new AllProductsCache(factors);
        // univariate multifactor diophantine solver
        UMultiDiophantineSolver<uPoly> uSolver = new UMultiDiophantineSolver<>(uFactors);

        UnivariatePolynomial<uPoly>[] solution = new UnivariatePolynomial[factors.length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = UnivariatePolynomial.constant(baseSeries.domain, uFactors.factors[i]);
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
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void bivariateLiftNoLCCorrection0(Poly base, Poly[] factors,
                                      IEvaluation<Term, Poly> evaluation, int degreeBound) {
        Domain<uPoly> uDomain = new UnivariatePolynomials<>((uPoly) factors[0].asUnivariate());
        UnivariatePolynomial<uPoly>[] res =
                bivariateLiftDense(seriesExpansionDense(uDomain, base, 1, evaluation),
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
            Term extends DegreeVector<Term>,
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
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void bivariateLift0(Poly base, Poly[] factors, Poly[] factorsLC,
                        IEvaluation<Term, Poly> evaluation, int degreeBound) {
        imposeLeadingCoefficients(factors, factorsLC);

        AllProductsCache<uPoly> uFactors = new AllProductsCache(asUnivariate(factors, evaluation));
        // univariate multifactor diophantine solver
        UMultiDiophantineSolver<uPoly> uSolver = new UMultiDiophantineSolver<>(uFactors);

        Domain<uPoly> uDomain = new UnivariatePolynomials<>(uFactors.get(0));
        UnivariatePolynomial<uPoly> baseSeries = seriesExpansionDense(uDomain, base, 1, evaluation);
        UnivariatePolynomial<uPoly>[] solution = new UnivariatePolynomial[factors.length];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = seriesExpansionDense(uDomain, factors[i], 1, evaluation);
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
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    uPoly[] asUnivariate(Poly[] array, IEvaluation<Term, Poly> evaluation) {
        uPoly u0 = (uPoly) evaluation.evaluateFrom(array[0], 1).asUnivariate();
        uPoly[] res = u0.arrayNewInstance(array.length);
        res[0] = u0;
        for (int i = 1; i < array.length; i++)
            res[i] = (uPoly) evaluation.evaluateFrom(array[i], 1).asUnivariate();
        return res;
    }

    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    Poly[] asMultivariate(uPoly[] array, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        Poly u0 = AMultivariatePolynomial.asMultivariate(array[0], nVariables, variable, ordering);
        Poly[] res = u0.arrayNewInstance(array.length);
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
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomial<uPoly> seriesExpansionDense(Domain<uPoly> domain, Poly poly, int variable, IEvaluation<Term, Poly> evaluate) {
        int degree = poly.degree(variable);
        uPoly[] coefficients = domain.createArray(degree + 1);
        for (int i = 0; i <= degree; i++)
            coefficients[i] = (uPoly) evaluate.taylorCoefficient(poly, variable, i).asUnivariate();
        return UnivariatePolynomial.create(domain, coefficients);
    }

    /**
     * Converts power series {@link #seriesExpansionDense(Domain, AMultivariatePolynomial, int, IEvaluation)}
     * back to multivariate polynomial
     */
    static <Term extends DegreeVector<Term>,
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
            Term extends DegreeVector<Term>,
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

            base = base.clone().multiply(CommonPolynomialsArithmetics.polyPow(lc, factors.length - 1, true));

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
            Term extends DegreeVector<Term>,
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
                factor.multiply(MultivariateReduction.divideExact(lcCorrection, factor.lc(0)));
            }

            base = base.clone().multiply(CommonPolynomialsArithmetics.polyPow(lc, factors.length - 1, true));

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
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    void multivariateLift0(Poly base,
                           Poly[] factors, Poly[] factorsLC,
                           IEvaluation<Term, Poly> evaluation,
                           int[] degreeBounds) {
        multivariateLift0(base, factors, factorsLC, evaluation, degreeBounds, 1);
    }

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
    static <Term extends DegreeVector<Term>,
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

        Domain<Poly> mDomain = new MultivariatePolynomials<>(factors[0]);
        for (int liftingVariable = from; liftingVariable < base.nVariables; ++liftingVariable) {
            Poly baseImage = evaluation.evaluateFrom(base, liftingVariable + 1);
            UnivariatePolynomial<Poly> baseSeries = seriesExpansion(mDomain, baseImage, liftingVariable, evaluation);
            UnivariatePolynomial<Poly>[] solution = new UnivariatePolynomial[factors.length];
            for (int i = 0; i < solution.length; i++) {
                solution[i] = seriesExpansion(mDomain, factors[i], liftingVariable, evaluation);
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

            if (liftingVariable < base.nVariables) {// don't perform on the last step
                dSolver.images[liftingVariable] = new AllProductsCache<>(
                        evaluation.evaluateFrom(factors, liftingVariable + 1))
                        .exceptArray();
            }
        }
    }

    /**
     * Generates a power series expansion for poly about the point specified by variable and evaluation
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    UnivariatePolynomial<Poly> seriesExpansion(Domain<Poly> domain, Poly poly, int variable, IEvaluation<Term, Poly> evaluate) {
        int degree = poly.degree(variable);
        Poly[] coefficients = domain.createArray(degree + 1);
        for (int i = 0; i <= degree; i++)
            coefficients[i] = evaluate.taylorCoefficient(poly, variable, i);
        return UnivariatePolynomial.create(domain, coefficients);
    }

    /**
     * Converts power series {@link #seriesExpansion(Domain, AMultivariatePolynomial, int, IEvaluation)}
     * back to multivariate polynomial
     */
    static <Term extends DegreeVector<Term>,
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

    static <Poly extends IGeneralPolynomial<Poly>>
    BernardinsTrick<Poly> createBernardinsTrick(UnivariatePolynomial<Poly>[] factors, int degreeBound) {
        if (Arrays.stream(factors).allMatch(UnivariatePolynomial::isConstant))
            return new BernardinsTrickWithoutLCCorrection<>(factors);
        else
            return new BernardinsTrickWithLCCorrection<>(factors, degreeBound);
    }

    static abstract class BernardinsTrick<Poly extends IGeneralPolynomial<Poly>> {
        // the factors
        final UnivariatePolynomial<Poly>[] factors;
        // partial products: {factor[0] * factor[1], (factor[0] * factor[1]) * factor[2], ... }
        final UnivariatePolynomial<Poly>[] partialProducts;
        final Domain<Poly> domain;

        BernardinsTrick(UnivariatePolynomial<Poly>... factors) {
            this.factors = factors;
            this.partialProducts = factors[0].arrayNewInstance(factors.length - 1);
            this.domain = factors[0].domain;

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
    static final class BernardinsTrickWithoutLCCorrection<Poly extends IGeneralPolynomial<Poly>>
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
            Poly newElement = domain.getZero();
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
                newElement = domain.getZero();
                for (int i = 1; i <= (pDegree + 1); i++)
                    newElement.add(partialProducts[j - 1].get(i).clone().multiply(factors[j + 1].get(pDegree - i + 1)));
                partialProducts[j].set(pDegree + 1, newElement);
            }
        }
    }

    /** Bernardin's trick for fast f_0 * f_1 * ... * f_N computing (leading coefficients are took into account) */
    static final class BernardinsTrickWithLCCorrection<Poly extends IGeneralPolynomial<Poly>>
            extends BernardinsTrick<Poly> {
        final int degreeBound;

        public BernardinsTrickWithLCCorrection(UnivariatePolynomial<Poly>[] factors, int degreeBound) {
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
}


