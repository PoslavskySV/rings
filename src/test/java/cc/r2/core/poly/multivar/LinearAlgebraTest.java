package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.util.ArraysUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.Collectors;

import static cc.r2.core.poly.multivar.LinearAlgebra.rowEchelonForm;
import static cc.r2.core.poly.multivar.LinearAlgebra.solve;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class LinearAlgebraTest {
    @Test
    public void test1() throws Exception {
        BigInteger[][] system = {
                {BigInteger.valueOf(1), BigInteger.valueOf(2), BigInteger.valueOf(3)},
                {BigInteger.valueOf(2), BigInteger.valueOf(2), BigInteger.valueOf(3)},
                {BigInteger.valueOf(3), BigInteger.valueOf(2), BigInteger.valueOf(3)},
        };
        BigInteger[] rhs = {
                BigInteger.valueOf(1),
                BigInteger.valueOf(2),
                BigInteger.valueOf(3),
        };


        BigInteger[] r = solve(new ModularDomain(SmallPrimes.nextPrime(12324)), system, rhs);
        System.out.println(Arrays.toString(r));
    }

    @Test
    public void test2() throws Exception {
        long[][] lhs0 = {
                {14945154, 0, 0, 0, 0, 0, 0, 0},
                {15840518, 0, 0, 0, 0, 0, 0, 23072526},
                {0, 1, 0, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0, 0, 18556231},
                {0, 0, 7274769, 0, 0, 0, 0, 0},
                {0, 0, 1285016, 0, 0, 0, 0, 2411651},
                {0, 0, 0, 9614891, 7274769, 0, 0, 0},
                {0, 0, 0, 4514307, 1285016, 0, 0, 17488741},
                {0, 0, 0, 0, 0, 9614891, 0, 0},
                {0, 0, 0, 0, 0, 4514307, 0, 7852752},
                {0, 0, 0, 0, 0, 0, 9614891, 0},
                {0, 0, 0, 0, 0, 0, 4514307, 22089485}};
        long[] rhs0 = {1, 0, 22655879, 0, 12777324, 0, 1128298, 0, 20152010, 0, 4506067, 0};
        long modulus = 23072527;

        BigInteger[][] lhs = convert(lhs0);
        BigInteger[] rhs = convert(rhs0);

        ModularDomain domain = new ModularDomain(modulus);
        BigInteger[] solution = solve(domain, lhs, rhs);
        long[] expected = {16402965, 22655879, 11505290, 2916536, 13894224, 7600529, 2132874, 14945154};
        Assert.assertArrayEquals(convert(expected), solution);
    }

    @Test
    public void test3() throws Exception {
        long[][] lhs0 = {
                {1, 0, 0, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 5642358},
                {0, 1, 1168354, 5331039, 0, 0, 0, 0},
                {0, 1, 798805, 1341857, 0, 0, 0, 5298367},
                {0, 0, 0, 0, 3103458, 0, 0, 0}, {0, 0, 0, 0, 1168354, 0, 0, 4274594},
                {0, 0, 0, 0, 0, 1, 1168354, 0}, {0, 0, 0, 0, 0, 1, 798805, 4627257}
        };
        long[] rhs0 = {1, 0, 740880, 0, 1693671, 0, 810986, 0};
        long modulus = 5642359;

        BigInteger[][] lhs = convert(lhs0);
        BigInteger[] rhs = convert(rhs0);

        ModularDomain domain = new ModularDomain(modulus);
        rowEchelonForm(domain, lhs, rhs);

        System.out.println(prettyMatrix(lhs0));
        System.out.println(prettyMatrix(lhs));
        System.out.println(Arrays.toString(rhs));

        BigInteger[] solution = solve(domain, lhs, rhs);
        long[] expected = {1, 2183072, 130178, 0, 6, 3367849, 8000, 1};
        Assert.assertArrayEquals(convert(expected), solution);
    }

    static BigInteger[] convert(long[] arr) {
        BigInteger[] r = new BigInteger[arr.length];
        for (int i = 0; i < arr.length; i++)
            r[i] = BigInteger.valueOf(arr[i]);
        return r;
    }

    static BigInteger[][] convert(long[][] arr) {
        BigInteger[][] r = new BigInteger[arr.length][];
        for (int i = 0; i < arr.length; i++)
            r[i] = convert(arr[i]);
        return r;
    }

    private static String padding(char c, int len) {
        return new String(ArraysUtil.arrayOf(c, len));
    }

    private static String padd(String str, int newLen) {
        return padding(' ', newLen - str.length()) + str;
    }

    public static String prettyMatrix(Object[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = matrix[i][j].toString();
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        return sb.append("}").toString();
    }

    public static String prettyMatrix(long[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = Long.toString(matrix[i][j]);
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        return sb.append("}").toString();
    }
}