package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.linear.LinearSolver;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.util.ArraysUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.MultivariateRing;
import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.GroebnerBasis;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;

/**
 * Utility methods based on Groebner bases
 */
public final class GroebnerMethods {
    private GroebnerMethods() {}

    /* *********************************************** Elimination *********************************************** */

    /**
     * Eliminates specified variables from the given ideal.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> eliminate(List<Poly> ideal, int... variables) {
        if (ideal.isEmpty())
            return Collections.emptyList();

        Comparator<DegreeVector> originalOrder = ideal.get(0).ordering;
        Comparator<DegreeVector> optimalOrder = optimalOrder(ideal);

        List<Poly> eliminationIdeal = ideal;
        if (!(optimalOrder instanceof MonomialOrder.GrevLexWithPermutation))
            for (int variable : variables)
                eliminationIdeal = GroebnerBasis(eliminationIdeal,
                        new MonomialOrder.EliminationOrder(optimalOrder, variable))
                        .stream()
                        .filter(p -> p.degree(variable) == 0)
                        .collect(Collectors.toList());
        else {
            MonomialOrder.GrevLexWithPermutation order = (MonomialOrder.GrevLexWithPermutation) optimalOrder;
            int[] inversePermutation = MultivariateGCD.inversePermutation(order.permutation);
            for (int variable : variables)
                eliminationIdeal = GroebnerBasis(eliminationIdeal
                        .stream()
                        .map(p -> AMultivariatePolynomial.renameVariables(p, order.permutation))
                        .collect(Collectors.toList()), new MonomialOrder.EliminationOrder(GREVLEX, inversePermutation[variable]))
                        .stream()
                        .map(p -> AMultivariatePolynomial.renameVariables(p, inversePermutation))
                        .filter(p -> p.degree(variable) == 0)
                        .collect(Collectors.toList());
        }

        return eliminationIdeal.stream().map(p -> p.setOrdering(originalOrder)).collect(Collectors.toList());
    }

    /* ******************************************* Algebraic dependence ******************************************** */

    /**
     * Returns true if a given set of polynomials is probably algebraically dependent or false otherwise (which means
     * that the given set is certainly independent). The method applies two criteria: it tests for lead set (LEX)
     * independence and does a probabilistic Jacobian test.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean probablyAlgebraicallyDependentQ(List<Poly> sys) {
        if (sys.isEmpty())
            return false;

        Poly representative = sys.get(0);
        if (sys.size() > representative.nVariables)
            return true;

        // give a check for LEX leading terms set
        List<DegreeVector> leadTerms;
        if (sys.stream().allMatch(p -> p.ordering == LEX))
            leadTerms = sys.stream().map(AMultivariatePolynomial::lt).collect(Collectors.toList());
        else
            leadTerms = sys.stream().map(p -> p.lt(LEX)).collect(Collectors.toList());
        if (!algebraicallyDependentMonomialsQ(leadTerms))
            return false;

        if (isMonomialIdeal(sys))
            return true;

        if (probablyMaximalJacobianRankQ(JacobianMatrix(sys)))
            return false;

        return true;
    }

    /**
     * Returns true if a given set of polynomials is algebraically dependent or false otherwise.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean algebraicallyDependentQ(List<Poly> sys) {
        return !algebraicRelations(sys).isEmpty();
    }

    /**
     * Tests for algebraic dependence the set of monomials
     */
    static boolean algebraicallyDependentMonomialsQ(List<DegreeVector> sys) {
        if (sys.isEmpty())
            return false;

        // build a homogeneous linear system
        int nVariables = sys.get(0).exponents.length;
        int nUnknowns = sys.size();

        // fixme use Bareiss in future
        @SuppressWarnings("unchecked")
        Rational<BigInteger>[][] lhs = new Rational[nVariables][nUnknowns];
        for (int i = 0; i < nVariables; ++i)
            for (int j = 0; j < nUnknowns; ++j)
                lhs[i][j] = Q.valueOf(sys.get(j).exponents[i]);
        Rational<BigInteger>[] rhs = Q.createZeroesArray(nVariables);

        // try to solve the system
        Rational<BigInteger>[] solution = Q.createZeroesArray(nUnknowns);
        LinearSolver.SystemInfo solveResult = LinearSolver.solve(Q, lhs, rhs, solution);
        if (solveResult == LinearSolver.SystemInfo.Consistent && Arrays.stream(solution).allMatch(Rational::isZero))
            return false;
        if (solveResult == LinearSolver.SystemInfo.Inconsistent)
            return false;
        return true;
    }

    /** Number of random substitutions for polynomial Jacobian to deduce its rank */
    private static final int N_JACOBIAN_EVALUATIONS_TRIES = 2;

    /** Probabilistic test for the maximality of the rank of Jacobian matrix */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean probablyMaximalJacobianRankQ(Poly[][] jacobian) {
        if (jacobian[0][0] instanceof MultivariatePolynomialZp64)
            return probablyMaximalJacobianRankQ((MultivariatePolynomialZp64[][]) jacobian);
        else
            return probablyMaximalJacobianRankQ((MultivariatePolynomial[][]) jacobian);
    }

    /** Probabilistic test for the maximality of the rank of Jacobian matrix */
    static boolean probablyMaximalJacobianRankQ(MultivariatePolynomialZp64[][] jacobian) {
        int nRows = jacobian.length, nColumns = jacobian[0].length;
        MultivariatePolynomialZp64 factory = jacobian[0][0];
        IntegersZp64 ring = factory.ring;
        long[][] matrix = new long[nRows][nColumns];
        long[] substitution = new long[nRows];
        RandomGenerator random = PrivateRandom.getRandom();
        for (int i = 0; i < N_JACOBIAN_EVALUATIONS_TRIES; ++i) {
            for (int var = 0; var < nRows; ++var)
                substitution[var] = ring.randomNonZeroElement(random);

            for (int iRow = 0; iRow < nRows; ++iRow)
                for (int iColumn = 0; iColumn < nColumns; ++iColumn)
                    matrix[iRow][iColumn] = jacobian[iRow][iColumn].evaluate(substitution);

            int nz = LinearSolver.rowEchelonForm(ring, matrix, null, false, true);
            if (nz == 0)
                return true;
        }
        return false;
    }

    /** Probabilistic test for the maximality of the rank of Jacobian matrix */
    static <E> boolean probablyMaximalJacobianRankQ(MultivariatePolynomial<E>[][] jacobian) {
        int nRows = jacobian.length, nColumns = jacobian[0].length;
        MultivariatePolynomial<E> factory = jacobian[0][0];
        Ring<E> ring = factory.ring;
        E[][] matrix = ring.createArray2d(nRows, nColumns);
        E[] substitution = ring.createArray(nRows);
        RandomGenerator random = PrivateRandom.getRandom();
        for (int i = 0; i < N_JACOBIAN_EVALUATIONS_TRIES; ++i) {
            for (int var = 0; var < nRows; ++var)
                substitution[var] = ring.randomNonZeroElement(random);

            for (int iRow = 0; iRow < nRows; ++iRow)
                for (int iColumn = 0; iColumn < nColumns; ++iColumn)
                    matrix[iRow][iColumn] = jacobian[iRow][iColumn].evaluate(substitution);

            // fixme use Bareiss in future
            int nz = LinearSolver.rowEchelonForm(ring, matrix, null, false, true);
            if (nz == 0)
                return true;
        }
        return false;
    }

    /**
     * Gives a list of algebraic relations (annihilating polynomials) for the given list of polynomials
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> algebraicRelations(List<Poly> polys) {
        if (!probablyAlgebraicallyDependentQ(polys))
            return Collections.emptyList();

        int nInitialVars = polys.get(0).nVariables;
        int nAdditionalVars = polys.size();
        List<Poly> helpPolys = new ArrayList<>();
        for (int i = 0; i < polys.size(); i++) {
            Poly p = polys.get(i).setNVariables(nInitialVars + nAdditionalVars);
            helpPolys.add(p.createMonomial(nInitialVars + i, 1).subtract(p));
        }

        int[] dropVars = ArraysUtil.sequence(0, nInitialVars);
        return eliminate(helpPolys, dropVars)
                .stream()
                .filter(p -> p.degree(dropVars) == 0)
                .map(p -> p.dropVariables(dropVars))
                .collect(Collectors.toList());
    }

    /**
     * Creates a Jacobian matrix of a given list of polynomials
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[][] JacobianMatrix(List<Poly> sys) {
        if (sys.isEmpty())
            throw new IllegalArgumentException("Empty list");
        MultivariateRing<Poly> ring = MultivariateRing(sys.get(0));
        Poly[][] jacobian = ring.createArray2d(ring.nVariables(), sys.size());
        for (int i = 0; i < ring.nVariables(); ++i)
            for (int j = 0; j < sys.size(); ++j)
                jacobian[i][j] = sys.get(j).derivative(i);
        return jacobian;
    }
}
