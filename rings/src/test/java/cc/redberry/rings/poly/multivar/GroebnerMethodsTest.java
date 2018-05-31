package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.MultivariateRing;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.multivar.GroebnerMethods.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static org.junit.Assert.*;

/**
 *
 */
public class GroebnerMethodsTest {

    @Test
    public void test1() {
        List<DegreeVector> dvs;

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}),
                new DegreeVector(new int[]{3, 4, 5}),
                new DegreeVector(new int[]{9, 8, 7}),
                new DegreeVector(new int[]{1, -1, -1}));
        assertTrue(algebraicallyDependentMonomialsQ(dvs));

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}));
        assertFalse(algebraicallyDependentMonomialsQ(dvs));

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}),
                new DegreeVector(new int[]{3, 1, 2}),
                new DegreeVector(new int[]{4, 1, 2}),
                new DegreeVector(new int[]{5, 1, 2}));
        assertTrue(algebraicallyDependentMonomialsQ(dvs));
    }


    @Test
    public void test2() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(3, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        assertFalse(probablyAlgebraicallyDependentQ(Arrays.asList(coder.parse("x"), coder.parse("y"))));
        assertFalse(probablyAlgebraicallyDependentQ(Arrays.asList(coder.parse("x^2 + x"), coder.parse("x^2 + y"))));
    }

    @Ignore
    @Test
    public void test3() {
        List<MultivariatePolynomial<Rational<BigInteger>>> katsura = GroebnerBasisData.katsura(8);
        eliminate(katsura.stream().map(p -> p.mapCoefficients(Zp64(17), r -> r.numerator().longValueExact())).collect(Collectors.toList()), 0)
                .stream().sorted(Comparator.comparingInt(AMultivariatePolynomial::nUsedVariables))
                .forEach(System.out::println);
    }

    @Test
    public void test4() {
        MultivariateRing<MultivariatePolynomialZp64> ring = Rings.MultivariateRingZp64(3, Zp64(32003));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "z");

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(
                coder.parse("x^2"),
                coder.parse("x*y"),
                coder.parse("y^5")
        );

        assertTrue(Ideal.create(eliminate(ideal, 0)).isPrincipal());
    }

    @Test
    public void test5() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(5, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "t", "s", "z");

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(
                coder.parse("x - t"),
                coder.parse("y - t^2"),
                coder.parse("z - t^3"),
                coder.parse("s - x + y^3")
        );

        Ideal eliminate = Ideal.create(eliminate(ideal, 2, 3), GREVLEX);
        Ideal el1 = Ideal.create(eliminate(eliminate(ideal, 2), 3), GREVLEX);
        Ideal el2 = Ideal.create(eliminate(eliminate(ideal, 3), 2), GREVLEX);

        assertEquals(eliminate, el1);
        assertEquals(eliminate, el2);
    }

    @Test
    public void test6() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(3, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "z");

        List<MultivariatePolynomial<Rational<BigInteger>>> sys;

        sys = Arrays.asList(
                coder.parse("x - y"),
                coder.parse("x + y"),
                coder.parse("x^2 + y^2"));
        assertNonEmptyAnnihilators(sys, algebraicRelations(sys));

        assertTrue(algebraicRelations(Arrays.asList(
                coder.parse("x * z"),
                coder.parse("y * z"))).isEmpty());

        sys = Arrays.asList(
                coder.parse("x * y + z^2"),
                coder.parse("z^2 + y^2"),
                coder.parse("x^2*y^2 -2 * x * y^3 + y^4"));
        assertNonEmptyAnnihilators(sys, algebraicRelations(sys));
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void assertAnnihilators(List<Poly> initialPolys, List<Poly> algebraicRelations) {
        assertTrue(algebraicRelations.stream().map(p -> p.composition(initialPolys)).allMatch(AMultivariatePolynomial::isZero));
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void assertNonEmptyAnnihilators(List<Poly> initialPolys, List<Poly> algebraicRelations) {
        assertFalse(algebraicRelations.isEmpty());
        assertAnnihilators(initialPolys, algebraicRelations);
    }
}