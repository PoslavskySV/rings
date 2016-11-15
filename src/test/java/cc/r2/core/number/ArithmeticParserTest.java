package cc.r2.core.number;

import org.junit.Test;

/**
 * Created by poslavsky on 01/11/2016.
 */
public class ArithmeticParserTest {
    @Test
    public void test1() throws Exception {
        ModBigIntegerRing ring = new ModBigIntegerRing(new BigInteger("12"));

//        ModBigInteger a = new ModBigInteger(ring, new BigInteger("14"));
//        ModBigInteger b = new ModBigInteger(ring, new BigInteger("2"));
//        ModBigInteger z = ring.getZero();
//        System.out.println(z);
//        System.out.println(a.add(b));

//        System.out.println(ArithmeticParser.parse("2*77", ring));
//        System.out.println(ArithmeticParser.parse("(1+2)*3", ring));
        System.out.println(ArithmeticParser.parse("1+(1+2)*(2*77+(1+2)*3-1)", ring));
    }
}