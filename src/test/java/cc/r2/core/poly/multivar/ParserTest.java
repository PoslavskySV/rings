package cc.r2.core.poly.multivar;

import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Rationals;
import gnu.trove.list.array.TIntArrayList;
import org.junit.Test;

import java.util.*;

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

    @Test
    public void test1() throws Exception {
        Rationals.Rationals.parse("+12");
    }

    @Test
    public void test2() throws Exception {
        System.out.println(Parser.parse("2/3*a*b^2 - 1/3*a^3*b^2", Rationals.Rationals, DegreeVector.LEX));
    }

}