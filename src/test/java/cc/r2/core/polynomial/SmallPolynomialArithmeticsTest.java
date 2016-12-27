package cc.r2.core.polynomial;

import org.junit.Test;

/**
 * Created by poslavsky on 27/12/2016.
 */
public class SmallPolynomialArithmeticsTest {
    @Test
    public void test1() throws Exception {

        MutableLongPoly a = MutableLongPoly.create(1, 4);
        MutableLongPoly b = MutableLongPoly.create(0, 2, 3);
        MutableLongPoly polyModulus = MutableLongPoly.create(0, 4, 0, 1);
        long modulus = 5;

        System.out.println(SmallPolynomialArithmetics.multiplyMod(a,b, polyModulus, modulus));

    }
}