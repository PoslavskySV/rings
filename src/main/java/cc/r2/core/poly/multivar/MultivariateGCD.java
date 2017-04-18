package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well1024a;

import java.util.HashSet;

import static cc.r2.core.poly.multivar.DivisionWithRemainderMultivariate.divideAndRemainder;
import static cc.r2.core.poly.multivar.DivisionWithRemainderMultivariate.dividesQ;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
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
    static MultivariatePolynomial<BigInteger> modularZpGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {
        return modularZpGCD(((ModularDomain) a.domain).modulus.longValue(), a, b, new RandomDataGenerator(new Well1024a()), a.nVariables - 1);
    }


    @SuppressWarnings("unchecked")
    static MultivariatePolynomial<BigInteger> modularZpGCD(
            long modulus,
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomDataGenerator rndd,
            int variable) {

        int nVariables = a.nVariables;
        if (variable == 0) {
            bMutablePolynomialZp poly = PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b));
            if (poly.degree() == 0)
                return a.createOne();
            MultivariatePolynomial<BigInteger> res = asMultivariate(poly, nVariables, variable, a.ordering);
            return res;
        }

        //content of a and b in Zp[x]
        MultivariatePolynomial<bMutablePolynomialZp>
                aConv = convertZp(a, variable),
                bConv = convertZp(b, variable);

        bMutablePolynomialZp
                aCont = PolynomialGCD(aConv.data.values()),
                bCont = PolynomialGCD(bConv.data.values());

        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(a, asMultivariate(aCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        a = qd[0];

        qd = divideAndRemainder(b, asMultivariate(bCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        b = qd[0];

        bMutablePolynomialZp
                c = PolynomialGCD(aCont, bCont),
                g = PolynomialGCD(aConv.lc(), bConv.lc());

        int n = Math.min(a.degrees()[variable - 1], b.degrees()[variable - 1]);

        Interpolation<BigInteger> interpolation = null;
        HashSet<BigInteger> seen = new HashSet<>();
        while (true) {
            if (seen.size() == modulus)
                throw new RuntimeException("too small field");
            BigInteger val = BigInteger.valueOf(rndd.nextLong(0, modulus - 1));
            if (seen.contains(val))
                continue;
            seen.add(val);
            BigInteger gVal = g.evaluate(val);
            if (gVal.isZero())
                continue;

            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, val),
                    bVal = b.evaluate(variable, val),
                    cVal = modularZpGCD(modulus, aVal, bVal, rndd, variable - 1);

            int m = cVal.degrees()[variable - 1];
            cVal = cVal.multiply(cVal.lc().modInverse(BigInteger.valueOf(modulus)).multiply(gVal));

            assert cVal.lc().equals(gVal);
            if (m == 0)
                return asMultivariate(c, nVariables, variable, a.ordering);

            if (m < n) {
                interpolation = new Interpolation<>(variable, val, cVal);
                n = m;
            } else if (m == n) {
                if (interpolation == null)
                    interpolation = new Interpolation<>(variable, val, cVal);
                else
                    interpolation.update(val, cVal);
            } else {
                int i = 0;
                continue;
            }

            if (interpolation.getPoints().size() <= 1)
                continue;

            MultivariatePolynomial<BigInteger> interpolated = fromZp(
                    convertZp(interpolation.getInterpolatingPolynomial(), variable).primitivePart(),
                    a.domain, variable);
            if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
                continue;

            return interpolated;
        }
    }
}
