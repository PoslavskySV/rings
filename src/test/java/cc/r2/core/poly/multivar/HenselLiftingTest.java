package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinatorialPort;
import cc.r2.core.combinatorics.IntTuplesPort;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.FiniteField;
import cc.r2.core.poly.MultivariatePolynomials;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.Homomorphism;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import cc.r2.core.util.TimeUnits;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
@Ignore
public class HenselLiftingTest {
    @Test
    public void teaaast1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(666));
        String[] vars = {"a", "b", "c", "d", "e"};

        Homomorphism homomorphism = new Homomorphism(new int[]{1, 2, 3, 4}, new long[]{0, 0, 0, 0});
        lMultivariatePolynomialZp
                root = parse("a^15 - 2*a*b^4 - 3*b*c + c + 2 + d^2*a - d*e", domain, vars),
                uSolution = root.evaluate(homomorphism.variables, homomorphism.values),
                square = root.clone().square();

        MultivariatePolynomials<lMultivariatePolynomialZp> mDomain = new MultivariatePolynomials<>(root.createOne());
        UnivariatePolynomial<lMultivariatePolynomialZp> equation =
                UnivariatePolynomial.create(mDomain, square, root.createZero(), root.createOne().negate());


        lMultivariatePolynomialZp lift = HenselLifting.lift(equation, homomorphism, uSolution);
        System.out.println(equation.evaluate(lift));
        System.out.println(lift);
        System.out.println(root);
    }

    @Test
    public void x() throws Exception {

        int[] set = new int[]{0, 1, 2};
        int nVariables = 4;
        int degree = 5;
        IntCombinatorialPort tups = new IntTuplesPort(ArraysUtil.arrayOf(nVariables, degree));
        int[] r;
        while ((r = tups.take()) != null)
            System.out.println(Arrays.toString(r));
    }

    @Test
    public void aatest2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(666));
        String[] vars = {"a", "b", "c", "d", "e"};
        Homomorphism homomorphism = new Homomorphism(new int[]{1, 2, 3, 4}, new long[]{0, 0, 0, 0});
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*b*c + c + 2 + d^2*a - d*e", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b*c + c + 2 + d*a - b*e", domain, vars),
                base = a.clone().multiply(b),
                ua = a.evaluate(homomorphism.variables, homomorphism.values),
                ub = b.evaluate(homomorphism.variables, homomorphism.values);

        System.out.println(MultivariateGCD.PolynomialGCD(a, b));

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            lMultivariatePolynomialZp[] r = HenselLifting.liftFactorization(base, homomorphism, ua, ub);
            Assert.assertEquals(a, r[0]);
            Assert.assertEquals(b, r[1]);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }


    }

    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                base = a.clone().multiply(b);

        System.out.println(MultivariateGCD.PolynomialGCD(a, b));

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftFactors0(base, aF, bF, aF.createOne(), bF.createOne());

//        System.out.println(a);
//        System.out.println(aF);
//
//        System.out.println(b);
//        System.out.println(bF);

        System.out.println(base.subtract(aF.clone().multiply(bF)));
    }

    @Test
    public void test1a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^2*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                base = a.clone().multiply(b);

        System.out.println(MultivariateGCD.PolynomialGCD(a, b));

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftFactors(base, aF, bF);

//        System.out.println(a);
//        System.out.println(aF);
//
//        System.out.println(b);
//        System.out.println(bF);

        System.out.println(base.subtract(aF.clone().multiply(bF)));
    }

    @Test
    public void test3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("1 + x - x^2 + x^2*y", domain, vars),
                b = parse("2 + x + x^2 + 2*y*x^2 + x^2*y^5 + y", domain, vars),
                base = a.clone().multiply(b);

        System.out.println(MultivariateGCD.PolynomialGCD(a, b));

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftFactors0(base,
                aF, bF,
                a.asUnivariate(0).lc(),
                b.asUnivariate(0).lc());
        System.out.println(base.clone().subtract(aF.clone().multiply(bF)));


        lMultivariatePolynomialZp baseLC = base.asUnivariate(0).lc();
        lMultivariatePolynomialZp nLC = baseLC.evaluate(1, 0);
        base = base.multiply(baseLC);
        aF = a.evaluate(1, 0).monic(nLC.lc());
        bF = b.evaluate(1, 0).monic(nLC.lc());

        HenselLifting.liftFactors0(base,
                aF, bF,
                baseLC,
                baseLC);

        System.out.println(base.subtract(aF.clone().multiply(bF)));
    }

    @Test
    public void test2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("1 + x + x^2 + x^2*y", domain, vars),
                b = parse("2 + x + 2*x^2 + x^2*y + y", domain, vars),
                base = a.clone().multiply(b);

        System.out.println(MultivariateGCD.PolynomialGCD(a, b));

        lMultivariatePolynomialZp lc = base.asUnivariate(0).lc();
        lMultivariatePolynomialZp a0 = a.clone(), b0 = b.clone(), base0 = base.clone();
//        a = a.multiply(lc);
//        b = b.multiply(lc);
//        base = base.multiply(lc.clone().square());

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftFactors(base, aF, bF);

        FiniteField<lUnivariatePolynomialZp> field = new FiniteField<>(lUnivariatePolynomialZ.monomial(1, base.degree(1)).modulus(domain));

        lUnivariatePolynomialZp
                aLC = aF.asUnivariate(0).lc().asUnivariate(),
                bLC = bF.asUnivariate(0).lc().asUnivariate();


        aF = aF.multiply(lMultivariatePolynomialZp.asMultivariate(field.reciprocal(aLC), base.nVariables, 1, base.ordering)).multiply(lc);
        for (lMonomialTerm term : aF.clone().terms) {
            if (term.exponents[1] >= base.degree(1))
                aF = aF.subtract(term);
        }
        lMultivariatePolynomialZp content = aF.asUnivariate(0).content();
        aF = AMultivariatePolynomial.asMultivariate(aF.asUnivariate(0).primitivePart(), 0).primitivePart();
        System.out.println(aF);
        System.out.println(a0);


        bF = bF.multiply(lMultivariatePolynomialZp.asMultivariate(field.reciprocal(bLC), base.nVariables, 1, base.ordering)).multiply(content);
        for (lMonomialTerm term : bF.clone().terms) {
            if (term.exponents[1] >= base.degree(1))
                bF = bF.subtract(term);
        }

        System.out.println(bF);
        System.out.println(b0);

        if (true) return;

        for (lMonomialTerm term : aF.clone().terms) {
            if (term.exponents[1] > base.degree(1))
                aF = aF.subtract(term);
        }

        for (lMonomialTerm term : bF.clone().terms) {
            if (term.exponents[1] > base.degree(1))
                bF = bF.subtract(term);
        }


        System.out.println(bF.asUnivariate(0).lc());
        System.out.println(lc);
//        assert (aF.asUnivariate(0).lc().equals(lc));
//        assert (bF.asUnivariate(0).lc().equals(lc));

        System.out.println(" ====  ");

        aF = AMultivariatePolynomial.asMultivariate(aF.asUnivariate(0).primitivePart(), 0).primitivePart();
        bF = AMultivariatePolynomial.asMultivariate(bF.asUnivariate(0).primitivePart(), 0).primitivePart();

        for (lMonomialTerm term : aF.clone().terms) {
            if (term.exponents[1] > base0.degree(1))
                aF = aF.subtract(term);
        }

        for (lMonomialTerm term : bF.clone().terms) {
            if (term.exponents[1] > base0.degree(1))
                bF = bF.subtract(term);
        }

        System.out.println(a0);
        System.out.println(aF);

        System.out.println(b0);
        System.out.println(bF);

//        System.out.println(a.toString(vars));
//        System.out.println(aF.toString(vars));
//
//        System.out.println(b.toString(vars));
//        System.out.println(bF.toString(vars));
//
//        System.out.println(base.subtract(aF.clone().multiply(bF)).toString(vars));
    }

    @Test
    public void asdasd() throws Exception {
        FiniteField<lUnivariatePolynomialZp> field = new FiniteField<>(lUnivariatePolynomialZ.create(0, 0, 0, 0, 1).modulus(7));

        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(1, 2, 1).modulus(7);
        lUnivariatePolynomialZp reciprocal = field.reciprocal(poly);

//        System.out.println(reciprocal);
        System.out.println(poly.clone().multiply(reciprocal)); ;
        System.out.println(field.multiply(poly, reciprocal));
//        System.out.println(reciprocal);

    }
}