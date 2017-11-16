package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariateSquareFreeFactorization;

import java.util.Arrays;

import static cc.redberry.rings.poly.multivar.Conversions64bit.asOverZp64;
import static cc.redberry.rings.poly.multivar.Conversions64bit.canConvertToZp64;
import static cc.redberry.rings.poly.multivar.MultivariateDivision.divideExact;

/**
 * @since 1.0
 */
public final class MultivariateSquareFreeFactorization {
    private MultivariateSquareFreeFactorization() {}

    /**
     * Performs square-free factorization of a {@code poly.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorization(Poly poly) {
        if (poly instanceof MultivariatePolynomial
                && ((MultivariatePolynomial) poly).ring.getZero() instanceof IPolynomial) {
            FactorDecomposition<Poly> factors = MultivariateFactorization.tryNested(poly,
                    MultivariateSquareFreeFactorization::SquareFreeFactorization);
            if (factors != null)
                return factors;
        }
        if (poly.coefficientRingCharacteristic().isZero())
            return SquareFreeFactorizationYunZeroCharacteristics(poly);
        else
            return SquareFreeFactorizationMusser(poly);
    }

    /**
     * Tests whether the given {@code poly} is square free.
     *
     * @param poly the polynomial
     * @return whether the given {@code poly} is square free
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean isSquareFree(Poly poly) {
        return MultivariateGCD.PolynomialGCD(poly, poly.derivative()).isConstant();
    }

    /**
     * Returns square-free part of the {@code poly}
     *
     * @param poly the polynomial
     * @return square free part
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly SquareFreePart(Poly poly) {
        return SquareFreeFactorization(poly).factors.stream().filter(x -> !x.isMonomial()).reduce(poly.createOne(), Poly::multiply);
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void addMonomial(FactorDecomposition<Poly> decomposition, Poly poly) {
        assert poly.isMonomial();
        decomposition.addConstantFactor(poly.lcAsPoly());
        poly = poly.monic();

        Term term = poly.lt();
        for (int i = 0; i < poly.nVariables; i++) {
            if (term.exponents[i] > 0)
                decomposition.addFactor(poly.create(term.singleVar(i)), term.exponents[i]);

        }
    }

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorUnivariate(Poly poly) {
        int var = poly.univariateVariable();
        return UnivariateSquareFreeFactorization
                .SquareFreeFactorization(poly.asUnivariate())
                .map(p -> AMultivariatePolynomial.asMultivariate((IUnivariatePolynomial) p, poly.nVariables, var, poly.ordering));
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] reduceContent(Poly poly) {
        Term monomialContent = poly.monomialContent();
        poly = poly.divideOrNull(monomialContent);

        Poly constantContent = poly.contentAsPoly();
        if (poly.signumOfLC() < 0)
            constantContent = constantContent.negate();

        poly = poly.divideByLC(constantContent);
        return poly.createArray(constantContent, poly.create(monomialContent));
    }

    /**
     * Performs square-free factorization of a {@code poly} which coefficient ring has zero characteristic using Yun's
     * algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationYunZeroCharacteristics(Poly poly) {
        if (!poly.coefficientRingCharacteristic().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        FactorDecomposition<Poly> decomposition
                = FactorDecomposition
                .constantFactor(content[0]);
        addMonomial(decomposition, content[1]);
        SquareFreeFactorizationYun0(poly, decomposition);
        return decomposition;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void SquareFreeFactorizationYun0(Poly poly, FactorDecomposition<Poly> factorization) {
        Poly[] derivative = poly.derivative();
        Poly gcd = MultivariateGCD.PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        Poly quot = divideExact(poly, gcd); // safely destroy (cloned) poly (not used further)
        Poly[] dQuot = poly.createArray(derivative.length);
        for (int i = 0; i < derivative.length; i++)
            dQuot[i] = divideExact(derivative[i], gcd);

        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            Poly[] qd = quot.derivative();
            for (int j = 0; j < derivative.length; j++)
                dQuot[j] = dQuot[j].subtract(qd[j]);

            Poly factor = MultivariateGCD.PolynomialGCD(quot, dQuot);
            quot = divideExact(quot, factor); // can destroy quot in divideAndRemainder
            for (int j = 0; j < derivative.length; j++)
                dQuot[j] = divideExact(dQuot[j], factor); // can destroy dQuot in divideAndRemainder
            if (!factor.isOne())
                factorization.addFactor(factor, i);
        }
    }

    /**
     * Performs square-free factorization of a {@code poly} which coefficient ring has zero characteristic using
     * Musser's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationMusserZeroCharacteristics(Poly poly) {
        if (!poly.coefficientRingCharacteristic().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        FactorDecomposition<Poly> decomposition
                = FactorDecomposition
                .constantFactor(content[0]);
        addMonomial(decomposition, content[1]);
        SquareFreeFactorizationMusserZeroCharacteristics0(poly, decomposition);
        return decomposition;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void SquareFreeFactorizationMusserZeroCharacteristics0(Poly poly, FactorDecomposition<Poly> factorization) {
        Poly[] derivative = poly.derivative();
        Poly gcd = MultivariateGCD.PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        Poly quot = divideExact(poly, gcd); // safely destroy (cloned) poly
        int i = 0;
        while (true) {
            ++i;
            Poly nextQuot = MultivariateGCD.PolynomialGCD(gcd, quot);
            gcd = divideExact(gcd, nextQuot); // safely destroy gcd (reassigned)
            Poly factor = divideExact(quot, nextQuot); // safely destroy quot (reassigned further)
            if (!factor.isConstant())
                factorization.addFactor(factor, i);
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
    }


    /**
     * Performs square-free factorization of a {@code poly} which coefficient ring has any characteristic using Musser's
     * algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationMusser(Poly poly) {
        if (poly.coefficientRingCharacteristic().isZero())
            throw new IllegalArgumentException("Positive characteristic expected");

        if (canConvertToZp64(poly))
            return SquareFreeFactorizationMusser(asOverZp64(poly)).map(Conversions64bit::convert);

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        Poly lc = poly.lcAsPoly();
        FactorDecomposition<Poly> fct = SquareFreeFactorizationMusser0(poly);
        addMonomial(fct, content[1]);
        return fct
                .addFactor(content[0], 1)
                .addFactor(lc, 1);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationMusser0(Poly poly) {
        poly.monic();
        if (poly.isConstant())
            return FactorDecomposition.constantFactor(poly);

        if (poly.degree() <= 1)
            return FactorDecomposition.singleFactor(poly.createOne(), poly);

        Poly[] derivative = poly.derivative();
        if (!Arrays.stream(derivative).allMatch(IPolynomial::isZero)) {
            Poly gcd = MultivariateGCD.PolynomialGCD(poly, derivative);
            if (gcd.isConstant())
                return FactorDecomposition.singleFactor(poly.createOne(), poly);
            Poly quot = divideExact(poly, gcd); // can safely destroy poly (not used further)

            FactorDecomposition<Poly> result = FactorDecomposition.constantFactor(poly.createOne());
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                Poly nextQuot = MultivariateGCD.PolynomialGCD(gcd, quot);
                Poly factor = divideExact(quot, nextQuot); // can safely destroy quot (reassigned further)
                if (!factor.isConstant())
                    result.addFactor(factor.monic(), i);

                gcd = divideExact(gcd, nextQuot); // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd);
                FactorDecomposition<Poly> gcdFactorization = SquareFreeFactorizationMusser0(gcd);
                gcdFactorization.raiseExponents(poly.coefficientRingCharacteristic().intValueExact());
                result.addAll(gcdFactorization);
                return result;
            } else
                return result;
        } else {
            Poly pRoot = pRoot(poly);
            FactorDecomposition<Poly> fd = SquareFreeFactorizationMusser0(pRoot);
            fd.raiseExponents(poly.coefficientRingCharacteristic().intValueExact());
            return fd.setConstantFactor(poly.createOne());
        }
    }

    /** p-th root of poly */
    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly pRoot(Poly poly) {
        if (poly instanceof MultivariatePolynomial)
            return (Poly) pRoot((MultivariatePolynomial) poly);
        else if (poly instanceof MultivariatePolynomialZp64)
            return (Poly) pRoot((MultivariatePolynomialZp64) poly);
        else
            throw new RuntimeException();
//        int modulus = poly.coefficientRingCharacteristic().intValueExact();
//        MonomialSet<Term> pRoot = new MonomialSet<>(poly.ordering);
//        for (Term term : poly) {
//            int[] exponents = term.exponents.clone();
//            for (int i = 0; i < exponents.length; i++) {
//                assert exponents[i] % modulus == 0;
//                exponents[i] = exponents[i] / modulus;
//            }
//            poly.add(pRoot, term.setDegreeVector(exponents));
//        }
//        return poly.create(pRoot);
    }

    private static <E> MultivariatePolynomial<E> pRoot(MultivariatePolynomial<E> poly) {
        Ring<E> ring = poly.ring;
        // p^(m -1) used for computing p-th root of elements
        BigInteger inverseFactor = ring.cardinality().divide(ring.characteristic());
        int modulus = poly.coefficientRingCharacteristic().intValueExact();
        MonomialSet<Monomial<E>> pRoot = new MonomialSet<>(poly.ordering);

        for (Monomial<E> term : poly) {
            int[] exponents = term.exponents.clone();
            for (int i = 0; i < exponents.length; i++) {
                assert exponents[i] % modulus == 0;
                exponents[i] = exponents[i] / modulus;
            }
            poly.add(pRoot, new Monomial<>(exponents, ring.pow(term.coefficient, inverseFactor)));
        }
        return poly.create(pRoot);
    }

    private static MultivariatePolynomialZp64 pRoot(MultivariatePolynomialZp64 poly) {
        assert !poly.ring.isPerfectPower();
        int modulus = MachineArithmetic.safeToInt(poly.ring.modulus);
        MonomialSet<MonomialZp64> pRoot = new MonomialSet<>(poly.ordering);

        for (MonomialZp64 term : poly) {
            int[] exponents = term.exponents.clone();
            for (int i = 0; i < exponents.length; i++) {
                assert exponents[i] % modulus == 0;
                exponents[i] = exponents[i] / modulus;
            }
            poly.add(pRoot, term.setDegreeVector(exponents));
        }
        return poly.create(pRoot);
    }
}
