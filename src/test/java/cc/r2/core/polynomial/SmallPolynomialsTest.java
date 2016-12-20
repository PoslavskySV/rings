package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static cc.r2.core.polynomial.SmallPolynomials.*;
import static org.junit.Assert.assertEquals;

public class SmallPolynomialsTest {
    @Test
    public void test1() throws Exception {
        MutableLongPoly a = new MutableLongPoly(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950);
        MutableLongPoly b = new MutableLongPoly(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216);

        SmallPolynomials.PolynomialRemainders prs = SmallPolynomials.Euclid(a, b, 11);
        System.out.println(prs.gcd());
    }

    @Test(expected = AssertionError.class)
    public void test2() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(new long[]{28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575, 0}, 7);
        MutableLongPoly divider = new MutableLongPoly(new long[]{24487310, 38204421, 12930314, 41553770, -1216266, 7382581, 15631547, 0, 0}, 6);
        pseudoDivideAndRemainder(dividend, divider);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(0, 14, 50, 93, 108, 130, 70);
        MutableLongPoly divider = new MutableLongPoly(63, 92, 143, 245, 146, 120, 90);
        MutableLongPoly gcd = PolynomialEuclid(dividend, divider, true).gcd();
        MutableLongPoly expected = new MutableLongPoly(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test4() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(7, -7, 0, 1);
        MutableLongPoly divider = new MutableLongPoly(-7, 0, 3);
        SmallPolynomials.PolynomialRemainders naive = PolynomialEuclid(dividend, divider, false);
        List<MutableLongPoly> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(new MutableLongPoly(63, -42));
        expectedNaive.add(new MutableLongPoly(-1L));
        assertEquals(expectedNaive, naive.remainders);

        SmallPolynomials.PolynomialRemainders primitive = PolynomialEuclid(dividend, divider, true);
        List<MutableLongPoly> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(new MutableLongPoly(-3, 2));
        expectedPrimitive.add(new MutableLongPoly(-1L));
        assertEquals(expectedPrimitive, primitive.remainders);


        SmallPolynomials.PolynomialRemainders subresultant = SubresultantEuclid(dividend, divider);
        List<MutableLongPoly> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(new MutableLongPoly(63, -42));
        expectedSubresultant.add(new MutableLongPoly(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }


    @Test
    public void test5() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(7, -7, 0, 1);
        MutableLongPoly divider = new MutableLongPoly(-7, 0, 3);
        SmallPolynomials.PolynomialRemainders a = PolynomialEuclid(dividend, divider, false);
        SmallPolynomials.PolynomialRemainders b = PolynomialEuclid(dividend, divider, true);
//        SmallPolynomials.PolynomialRemainders c = EuclideanPRS(dividend, divider);

        System.out.println(a.remainders);
        System.out.println(b.remainders);
        System.out.println(SubresultantEuclid(dividend, divider).remainders);

    }


    private static final int CBOUND = 1000;

    public static MutableLongPoly randomPoly(int degree, RandomGenerator rnd) {
        long[] data = new long[degree + 1];
        for (int i = 0; i <= degree; ++i) {
            data[i] = rnd.nextInt(CBOUND);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = -data[i];
        }
        while (data[degree] == 0)
            data[degree] = rnd.nextInt(CBOUND);
        return new MutableLongPoly(data, degree);
    }
}