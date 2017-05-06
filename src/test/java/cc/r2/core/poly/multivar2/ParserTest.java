package cc.r2.core.poly.multivar2;

import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class ParserTest {
    @Test
    public void name() throws Exception {
//        TIntArrayList ints = new TIntArrayList(new int[]{0, 1, 2, 3, 4, 5, 6});
//        ints.removeA(0);
//        System.out.println(ints);
        System.out.println(Parser.parse("a^2 * b^33 + a * c", MonomialTerm.GREVLEX));
    }
}