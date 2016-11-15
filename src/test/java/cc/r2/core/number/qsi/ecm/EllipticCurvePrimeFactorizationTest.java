package cc.r2.core.number.qsi.ecm;

import cc.r2.core.number.BigInteger;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by poslavsky on 14/11/2016.
 */
public class EllipticCurvePrimeFactorizationTest {
    @Test
    public void name1() throws Exception {
//        int[] queue = new int[200];
//        int squfof = Siqs.SQUFOF(1222222222L, queue);
//        for (int i = 0; i < 100; i++) {
//            System.out.println(squfof);
//        }
//
//        System.out.println(Arrays.toString(queue));
        ecm els = new ecm();
        els.NumberToFactor = new BigInteger("125");
        els.init();
        els.startNewFactorization(false);
//
        System.out.println(els.NbrFactors);
        System.out.println(Arrays.toString(els.PD));
        System.out.println(Arrays.toString(els.Exp1));

    }
}