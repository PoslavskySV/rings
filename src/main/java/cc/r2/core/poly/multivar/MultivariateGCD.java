package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.univar.PolynomialGCD.PolynomialGCD;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateGCD {
    private MultivariateGCD() {}

    private static void checkModGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {

        if (!(a.domain instanceof ModularDomain))
            throw new IllegalArgumentException();

    }

    @SuppressWarnings("unchecked")
    static MultivariatePolynomial<BigInteger> denseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {
        BigInteger domainSize = a.domain.size();
        if (domainSize.isInt() && Math.min(a.degree(), b.degree()) > domainSize.intValue())
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");

        int evaluationStackLimit = domainSize.isInt() ? domainSize.intValue() : -1;
        int[] aDegrees = a.degrees(),
                bDegrees = b.degrees(),
                degreeBounds = new int[a.nVariables];
        for (int i = 0; i < a.nVariables; i++)
            degreeBounds[i] = Math.min(aDegrees[i], bDegrees[i]);

        MultivariatePolynomial<BigInteger> result = denseModularGCD(a, b, PrivateRandom.getRandom(), a.nVariables - 1, degreeBounds, evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return result;
    }

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> denseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("input is zero");

        MultivariatePolynomial<BigInteger> factory = a;
        if (a.isConstant() || b.isConstant())
            return factory.createOne();

        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            bMutablePolynomialZp gcd = PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b));
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        Domain<BigInteger> domain = a.domain;

        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        MultivariatePolynomial<bMutablePolynomialZp>
                aConv = convertZp(a, variable),
                bConv = convertZp(b, variable);

        // content of a and b in Zp[x_k]
        bMutablePolynomialZp
                aCont = PolynomialGCD(aConv.data.values()),
                bCont = PolynomialGCD(bConv.data.values());

        // normalize a and b
        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(a, asMultivariate(aCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        a = qd[0];

        qd = divideAndRemainder(b, asMultivariate(bCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        b = qd[0];

        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = PolynomialGCD(aCont, bCont),
                lcGCD = PolynomialGCD(aConv.lc(), bConv.lc());

//        List<bMutablePolynomialZp> evalTest = new ArrayList<>();
//        for (int i = 0; i < variable; ++i)
//            evalTest.add(PolynomialGCD(convertZp(a, i).lc(), convertZp(b, i).lc()));

        int prevVarExponent = degreeBounds[variable - 1];
        Interpolation<BigInteger> interpolation = null;
        MultivariatePolynomial<BigInteger> previousInterpolation;
        Set<BigInteger> evaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // do division check (last chance)
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            BigInteger randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;

            evaluationStack.add(randomPoint);

            BigInteger lcVal = lcGCD.evaluate(randomPoint);
            if (lcVal.isZero())
                continue;

//            for (bMutablePolynomialZp lc : evalTest)
//                if (lc.evaluate(randomPoint).isZero())
//                    continue main;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<BigInteger> cVal = denseModularGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            // normalize gcd
            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            if (currExponent < prevVarExponent) {
                //better degree bound detected => start over
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
                continue;
            }

            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                continue;
            }

            previousInterpolation = interpolation.getInterpolatingPolynomial();
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] == 0
                    || (previousInterpolation != null && previousInterpolation.equals(interpolation.getInterpolatingPolynomial()))) {
                MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    private static MultivariatePolynomial<BigInteger> doDivisionCheck(
            MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, bMutablePolynomialZp contentGCD,
            Interpolation<BigInteger> interpolation, int variable) {
        if (interpolation == null)
            return null;
        MultivariatePolynomial<BigInteger> interpolated =
                fromZp(convertZp(interpolation.getInterpolatingPolynomial(), variable).primitivePart(), a.domain, variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        return interpolated.multiply(asMultivariate(contentGCD, a.nVariables, variable, a.ordering));
    }
}
