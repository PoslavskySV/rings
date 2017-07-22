package cc.r2.core.poly.multivar;

import cc.r2.core.poly.FiniteField;
import cc.r2.core.poly.Rationals;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

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

    @Test
    @SuppressWarnings("unchecked")
    public void test3() throws Exception {
        String[] strs = {
                "b+(1+x^2)*b^2+(x)*b^3+(1+x+x^2)*b^4",
                "(1+x)*b^4+(1+x)*b^8",
                "(x)*b^3+(1+x^2)*b^4+(x+x^2)*b^5",
                "(x+x^2)*b^2+(x^2)*b^3+(x+x^2)*b^4+(1+x^2)*b^5+(1+x^2)*b^6",
                "1+(1+x^2)*b+(x+x^2)*b^3+(x^2)*b^4+(x)*b^6+(1+x+x^2)*b^7",
                "(1+x+x^2)*b^7",
                "(1+x)*b^5",
                "(1+x)*b^4+(1+x)*b^8",
                "(1+x^2)*b^4+b^5+(x)*b^6+(1+x)*b^7",
                "(x)*b^3+(1+x^2)*b^4+(x+x^2)*b^5",
                "(x)*b^2+(1+x^2)*b^3+(1+x+x^2)*b^5+b^6+(1+x^2)*b^7+(x+x^2)*b^8",
                "(x+x^2)*b^2+(x^2)*b^3+(x+x^2)*b^4+(1+x^2)*b^5+(1+x^2)*b^6",
                "(1+x+x^2)*b^7",
                "1+(1+x^2)*b+(x+x^2)*b^3+(x^2)*b^4+(x)*b^6+(1+x+x^2)*b^7",
                "(1+x)+(x)*b+b^3+(1+x^2)*b^4",
                "(1+x+x^2)+b",
                "(x)*b^3+(1+x^2)*b^4+(x+x^2)*b^5",
                "(1+x+x^2)*b^7",
                "(1+x)*b^5",
                "(1+x)*b^4+(1+x)*b^8",
                "(1+x^2)*b^4+b^5+(x)*b^6+(1+x)*b^7",
                "(x)*b^3+(1+x^2)*b^4+(x+x^2)*b^5",
                "(x)*b^2+(1+x^2)*b^3+(1+x+x^2)*b^5+b^6+(1+x^2)*b^7+(x+x^2)*b^8",
                "(1+x+x^2)*b^7",
                "(1+x)*b^4+(1+x)*b^8",
        };
        FiniteField<lUnivariatePolynomialZp> minorDomain = new FiniteField<>(lUnivariatePolynomialZ.create(1, 0, 1, 1).modulus(2));
        FiniteField<UnivariatePolynomial<lUnivariatePolynomialZp>> domain = new FiniteField<>(UnivariatePolynomial.parse(minorDomain, "(1+x^2)+(x^2)*x+(x+x^2)*x^2+x^3"));
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<lUnivariatePolynomialZp> arr[] = Arrays.stream(strs)
                .map(s -> MultivariatePolynomial.parse(s, domain, vars))
                .toArray(MultivariatePolynomial[]::new);
        for (int i = 0; i < arr.length; i++) {
            Assert.assertEquals(strs[i], arr[i].toString());
            Assert.assertEquals(arr[i], arr[i].parsePoly(arr[i].toString()));
        }
    }
}