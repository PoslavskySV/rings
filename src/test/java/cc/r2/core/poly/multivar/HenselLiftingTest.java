package cc.r2.core.poly.multivar;

import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.Evaluation;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.HenselLifting.liftWang;
import static cc.r2.core.poly.multivar.HenselLifting.modImage;
import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
@Ignore
public class HenselLiftingTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^2*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("1 + x - x^2 + x^2*y", domain, vars),
                b = parse("2 + x + x^2 + 2*y*x^2 + x^2*y^5 + y", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 + a^15*b^2*c^3 + a^15*c + 7*a^15- 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5  + a^5*b*c - a^5*b + 2*a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();


        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test6() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                b = parse("1227874+3587706*b+5373508*a+7197578*a^2+a^3", domain, vars),
                a = parse("9540707+24*a", domain, vars),
                base = parse("7717493+597721*b+6517458*b^2+361611*a+9241048*a*b+9607947*a*b^2+3165308*a^2+9338813*a^3+24*a^4", domain, vars);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }


    @Test
    public void test8() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};

        int[] svars = {1, 2};
        long[] subs = {11, 12};
        Evaluation evaluation = new Evaluation(vars.length, subs, domain, DegreeVector.LEX);
        lMultivariatePolynomialZp
                a = parse("1 + c*b + b*a*c^5 + b*c*2*a^2", domain, vars).square().square(),
                b = parse("1227874+3587706*b+5373508*a+7197578*c^2*a^2+a^3", domain, vars).square().square(),
                base = a.clone().multiply(b),
                ua = a.evaluate(svars, subs),
                ub = b.evaluate(svars, subs);

        System.out.println(ua);
        System.out.println(a);
        HenselLifting.liftWang(base, ua, ub, evaluation);

        System.out.println(ua.subtract(a));
        System.out.println(ub.subtract(b));

    }

    @Test
    public void asd() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        lUnivariatePolynomialZp
                a = parse("1+11*a+24*a^2", domain).asUnivariate(),
                b = parse("2260692+5373508*a+7197578*a^2+a^3", domain).asUnivariate(),
                rhs = parse("1284104+2030309*a+8124550*a^2+7657516*a^3", domain).asUnivariate();

        lUnivariatePolynomialZp[] sol = HenselLifting.solveDiophantine(a, b, rhs);
        System.out.println(a.multiply(sol[0]).add(b.multiply(sol[1])).subtract(rhs));
//        System.out.println(UnivariateGCD.PolynomialGCD(a, b));
        //1284104+2030309*a+8124550*a^2+7657516*a^3

    }

    @Test
    public void testHenselLiftingRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = 1000;
        for (int i = 0; i < nIterations; i++) {
            System.out.println(i);
            long modulus = getModulusRandom(rndd.nextInt(5, 30));
            lIntegersModulo domain = new lIntegersModulo(modulus);

            int nVars = rndd.nextInt(2, 3);

            lMultivariatePolynomialZp[] polys = new lMultivariatePolynomialZp[2];
            for (int j = 0; j < polys.length; j++) {
                polys[j] = MultivariatePolynomial.asLongPolyZp(
                        RandomMultivariatePolynomial.randomPolynomial(nVars,
                                rndd.nextInt(5, 15), rndd.nextInt(5, 25), domain.asDomain(), DegreeVector.LEX, rnd));
            }

            lMultivariatePolynomialZp a = polys[0], b = polys[1];
//            //make monic
//            a = a.add(a.createUnivariateMonomial(0, a.degree(0) + 2));
//            b = b.add(b.createUnivariateMonomial(0, b.degree(0) + 2));

            System.out.println("gcd");
            if (!MultivariateGCD.PolynomialGCD(a, b).isOne()) {
                System.out.println("bad inp");
                --i;
                continue;
            }
            System.out.println("gcd...done");

            int nEvaluations = 10;
            evals:
            for (int n = 0; n < nEvaluations; n++) {
                long[] substitutions = new long[nVars - 1];
                for (int j = 0; j < substitutions.length; j++) {
                    do {
                        substitutions[j] = domain.randomElement(rnd);
                    } while (substitutions[j] == 0);
                }


                Evaluation evaluation = new Evaluation(nVars, substitutions, domain, DegreeVector.LEX);
                lMultivariatePolynomialZp
                        ua = evaluation.evaluateFrom(a, 1),
                        ub = evaluation.evaluateFrom(b, 1);
                if (!MultivariateGCD.PolynomialGCD(ua, ub).isOne()) {
                    System.out.println("bad p");
                    --n;
                    continue;
                }

                if (!ua.getSkeleton().equals(a.getSkeleton(0)) || !ub.getSkeleton().equals(b.getSkeleton(0))) {
                    System.out.println("bad p");
                    --n;
                    continue;
                }

                try {
                    //System.out.println(a);
                    //System.out.println(b);
                    liftWang(a.clone().multiply(b), ua, ub, a.lc(0), b.lc(0), evaluation);
                    Assert.assertEquals(a, ua);
                    Assert.assertEquals(b, ub);
                } catch (Throwable thr) {
                    System.out.println(domain);
                    System.out.println(a);
                    System.out.println(b);
                    System.out.println(Arrays.toString(evaluation.values));
                    throw thr;
                }
            }
        }
    }

    @Test
    public void testHenselLifting1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};

        lMultivariatePolynomialZp
                a = parse("1 + c*b + b*a*c^5 + b*c*2*a^2", domain, vars),
                b = parse("1227874+3587706*b+5373508*a+7197578*c^2*a^2+a^3", domain, vars),
                base = a.clone().multiply(b);

        long[] subs = {11, 12};
        Evaluation evaluation = new Evaluation(vars.length, subs, domain, DegreeVector.LEX);

        lMultivariatePolynomialZp
                ua = evaluation.evaluateFrom(a, 1),
                ub = evaluation.evaluateFrom(b, 1);

        HenselLifting.liftWang(base, ua, ub, a.lc(0), b.lc(0), evaluation);

        Assert.assertEquals(a, ua);
        Assert.assertEquals(b, ub);
    }

    @Test
    public void testHenselLifting2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(31012727);
        String[] vars = {"a", "b"};

        lMultivariatePolynomialZp
                a = parse("26590303+a^1", domain, vars),
                b = parse("24266411+16065264*a^1+8307958*a^2", domain, vars),
                base = parse("31012723*b^2+30996897*b^3+30965273*b^4+7876*a*b+71042*a*b^2+31009448*a*b^3+31002926*a*b^4+30989135*a^2+12*a^2*b+31008764*a^2*b^2+30971435*a^2*b^3+41292*a^3*b", domain, vars),
        sf1 = parse("-9801", domain, vars),
        sf2 = parse("b^2-a+10337576b", domain, vars),
        sf3 = parse("5496291*a^2*b+a*b^2-10840690*a*b-9919891*b^2+7053091*a-10840690*b", domain, vars);

        long[] subs = {3793049};
        Evaluation evaluation = new Evaluation(vars.length, subs, domain, DegreeVector.LEX);

        lMultivariatePolynomialZp
                ua = evaluation.evaluateFrom(a, 1),
                ub = evaluation.evaluateFrom(b, 1);

        HenselLifting.liftWang(base, ua, ub, null, null, evaluation);

        System.out.println(base.clone().evaluate(1, subs[0]).subtract(a.clone().multiply(b)));
        System.out.println(base.clone().subtract(ua.clone().multiply(ub)));
        System.out.println(base.clone().subtract(ua.clone().multiply(ub)).evaluate(1, subs[0]));
//        Assert.assertEquals(a, ua);
//        Assert.assertEquals(b, ub);
    }
}