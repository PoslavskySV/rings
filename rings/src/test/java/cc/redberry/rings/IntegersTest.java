package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import org.junit.Test;

import static cc.redberry.rings.Rings.Z;
import static org.junit.Assert.*;

/**
 * @author Stanislav Poslavsky
 * @since 2.2
 */
public class IntegersTest {
    @Test
    public void test1() throws Exception {
        BigInteger p = Z.valueOf(2 * 3 * 4 * 5);
        System.out.println(Z.factor(p));
    }
}