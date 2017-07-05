package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.poly.LongArithmetics;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.MultivariateReduction.divideExact;

/**
 * @author Stanislav Poslavsky
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
        if (poly.coefficientDomainCharacteristics().isZero())
            return SquareFreeFactorizationYunZeroCharacteristics(poly);
        else
            return SquareFreeFactorizationMusser(poly);
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] reduceContent(Poly poly) {
        Term monomialContent = poly.monomialContent();
        poly = poly.divideOrNull(monomialContent);

        Poly constantContent = poly.contentAsPoly();
        if (poly.signum() < 0)
            constantContent = constantContent.negate();

        poly = poly.divideByLC(constantContent);
        return poly.arrayNewInstance(constantContent, poly.create(monomialContent));
    }

    @SuppressWarnings("unchecked")
    private static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Poly poly, Poly[] arr) {
        Poly[] all = (Poly[]) poly.arrayNewInstance(arr.length + 1);
        all[0] = poly;
        System.arraycopy(arr, 0, all, 1, arr.length);
        return MultivariateGCD.PolynomialGCD(all);
    }

    /**
     * Performs square-free factorization of a {@code poly} which coefficient domain has zero characteristic
     * using Yun's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationYunZeroCharacteristics(Poly poly) {
        if (!poly.coefficientDomainCharacteristics().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        FactorDecomposition<Poly> decomposition
                = FactorDecomposition
                .constantFactor(content[0])
                .addFactor(content[1], 1);
        SquareFreeFactorizationYun0(poly, decomposition);
        return decomposition;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void SquareFreeFactorizationYun0(Poly poly, FactorDecomposition<Poly> factorization) {
        Poly[] derivative = poly.derivative();
        Poly gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        Poly quot = divideExact(poly, gcd); // safely destroy (cloned) poly (not used further)
        Poly[] dQuot = poly.arrayNewInstance(derivative.length);
        for (int i = 0; i < derivative.length; i++)
            dQuot[i] = divideExact(derivative[i], gcd);

        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            Poly[] qd = quot.derivative();
            for (int j = 0; j < derivative.length; j++)
                dQuot[j] = dQuot[j].subtract(qd[j]);

            Poly factor = PolynomialGCD(quot, dQuot);
            quot = divideExact(quot, factor); // can destroy quot in divideAndRemainder
            for (int j = 0; j < derivative.length; j++)
                dQuot[j] = divideExact(dQuot[j], factor); // can destroy dQuot in divideAndRemainder
            if (!factor.isOne())
                factorization.addFactor(factor, i);
        }
    }

    /**
     * Performs square-free factorization of a {@code poly} which coefficient domain has zero characteristic
     * using Musser's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationMusserZeroCharacteristics(Poly poly) {
        if (!poly.coefficientDomainCharacteristics().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        FactorDecomposition<Poly> decomposition
                = FactorDecomposition
                .constantFactor(content[0])
                .addFactor(content[1], 1);
        SquareFreeFactorizationMusserZeroCharacteristics0(poly, decomposition);
        return decomposition;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void SquareFreeFactorizationMusserZeroCharacteristics0(Poly poly, FactorDecomposition<Poly> factorization) {
        Poly[] derivative = poly.derivative();
        Poly gcd = PolynomialGCD(poly, derivative);
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
     * Performs square-free factorization of a {@code poly} which coefficient domain has non zero characteristic
     * using Musser's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> SquareFreeFactorizationMusser(Poly poly) {
        if (poly.coefficientDomainCharacteristics().isZero())
            throw new IllegalArgumentException("Positive characteristics expected");

        poly = poly.clone();
        Poly[] content = reduceContent(poly);
        Poly lc = poly.lcAsPoly();
        return SquareFreeFactorizationMusser0(poly)
                .addFactor(content[1], 1)
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
        if (!Arrays.stream(derivative).allMatch(IGeneralPolynomial::isZero)) {
            Poly gcd = PolynomialGCD(poly, derivative);
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
                gcdFactorization.raiseExponents(poly.coefficientDomainCharacteristics().intValueExact());
                result.addAll(gcdFactorization);
                return result;
            } else
                return result;
        } else {
            Poly pRoot = pRoot(poly);
            FactorDecomposition<Poly> fd = SquareFreeFactorizationMusser0(pRoot);
            fd.raiseExponents(poly.coefficientDomainCharacteristics().intValueExact());
            return fd.setConstantFactor(poly.createOne());
        }
    }

    /** p-th root of poly */
    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly pRoot(Poly poly) {
        if (poly instanceof MultivariatePolynomial)
            return (Poly) pRoot((MultivariatePolynomial) poly);
        else if (poly instanceof lMultivariatePolynomialZp)
            return (Poly) pRoot((lMultivariatePolynomialZp) poly);
        else
            throw new RuntimeException();
//        int modulus = poly.coefficientDomainCharacteristics().intValueExact();
//        MonomialsSet<Term> pRoot = new MonomialsSet<>(poly.ordering);
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
        Domain<E> domain = poly.domain;
        // p^(m -1) used for computing p-th root of elements
        BigInteger inverseFactor = domain.cardinality().divide(domain.characteristics());
        int modulus = poly.coefficientDomainCharacteristics().intValueExact();
        MonomialsSet<MonomialTerm<E>> pRoot = new MonomialsSet<>(poly.ordering);

        for (MonomialTerm<E> term : poly) {
            int[] exponents = term.exponents.clone();
            for (int i = 0; i < exponents.length; i++) {
                assert exponents[i] % modulus == 0;
                exponents[i] = exponents[i] / modulus;
            }
            poly.add(pRoot, new MonomialTerm<>(exponents, domain.pow(term.coefficient, inverseFactor)));
        }
        return poly.create(pRoot);
    }

    private static lMultivariatePolynomialZp pRoot(lMultivariatePolynomialZp poly) {
        assert !poly.domain.isPerfectPower();
        int modulus = LongArithmetics.safeToInt(poly.domain.modulus);
        MonomialsSet<lMonomialTerm> pRoot = new MonomialsSet<>(poly.ordering);

        for (lMonomialTerm term : poly) {
            int[] exponents = term.exponents.clone();
            for (int i = 0; i < exponents.length; i++) {
                assert exponents[i] % modulus == 0;
                exponents[i] = exponents[i] / modulus;
            }
            poly.add(pRoot, term.setDegreeVector(exponents));
        }
        return poly.create(pRoot);
    }

    /**
     * Tests whether the given {@code poly} is square free.
     *
     * @param poly the polynomial
     * @return whether the given {@code poly} is square free
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean isSquareFree(Poly poly) {
        return PolynomialGCD(poly, poly.derivative()).isConstant();
    }
}
