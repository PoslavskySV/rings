package cc.r2.core.polynomial;

import org.junit.Test;

import static cc.r2.core.polynomial.SmallPolynomialArithmetics.*;
import static org.junit.Assert.assertEquals;

public class SmallPolynomialArithmeticsTest {
    @Test
    public void test1() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1, 4);
        MutableLongPoly b = MutableLongPoly.create(0, 2, 3);
        MutableLongPoly polyModulus = MutableLongPoly.create(0, 4, 0, 1);
        long modulus = 5;
        assertEquals(MutableLongPoly.create(0, 4, 1), multiplyMod(a, b, polyModulus, modulus));
    }

    @Test
    public void test2() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1, 2, 3, 4, 0);
        assertEquals(MutableLongPoly.create(2, 6, 12), derivative(a));
    }

    @Test
    public void test3() throws Exception {
        assertEquals(MutableLongPoly.create(1, 2, 1), pow(MutableLongPoly.create(1, 1), 2));
    }

    @Test
    public void test4() throws Exception {
//        System.out.println(LongArithmetics.modInverse(7, 7));
//        assertEquals(MutableLongPoly.create(1, 2, 1), pow(MutableLongPoly.create(1, 1), 2));
        MutableLongPoly a = MutableLongPoly.create(1, 0, 1, 0, 1);
        MutableLongPoly b = MutableLongPoly.create(1, 1, 1);
        System.out.println(SmallPolynomialArithmetics.powMod(b,2,2));
    }
}