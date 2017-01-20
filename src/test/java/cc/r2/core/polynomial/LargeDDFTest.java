package cc.r2.core.polynomial;

import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.LongArithmetics.symMod;
import static cc.r2.core.polynomial.MutablePolynomial.createMonomial;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyMod;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyPowMod;
import static cc.r2.core.polynomial.DivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.DivideAndRemainder.quotient;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 09/01/2017.
 */
@Ignore
public class LargeDDFTest {
    static final MutablePolynomial bigPoly = MutablePolynomial.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1);
    static final int bigModulus = 5659;


    public static ArrayList<MutablePolynomial> xPowers(MutablePolynomial poly, InverseModMonomial invMod, long modulus) {
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        exps.add(MutablePolynomial.one());
        MutablePolynomial base = polyMod(createMonomial(1, (int) modulus), poly, invMod, modulus, false);
        exps.add(base);
        MutablePolynomial prev = base;
        for (int i = 0; i < poly.degree; i++)
            exps.add(prev = polyMod(prev.clone().multiply(base, modulus), poly, invMod, modulus, false));
        return exps;
    }

    @Test
    public void testXPowers1() throws Exception {
        long modulus = 7;
        MutablePolynomial poly = MutablePolynomial.create(6, 2, 3, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        ArrayList<MutablePolynomial> xPowers = xPowers(poly, invMod, modulus);

        for (int i = 0; i < xPowers.size(); i++) {
            System.out.println(xPowers.get(i));
        }
//        int m = dividend.degree - divider.degree;
//        MutablePolynomial q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1), modulus), m + 1, false).reverse();
//        if (q.degree < m)
//            q.shiftRight(m - q.degree);

        System.out.println("------");
//        for (int i = 0; i < xPowers.size(); i++) {
//            int monomialExponent  = i * (int) modulus;
//            int m = monomialExponent  - poly.degree;
//            if (m > 0) {
//                MutablePolynomial q = remainderMonomial(invMod.getInverse(m + 1), m + 1, false).reverse();
//                if (q.degree < m)
//                    q.shiftRight(m - q.degree);
//                q.multiply(poly, modulus).negate(modulus);
//                q.ensureCapacity(monomialExponent);
//                q.data[monomialExponent] += 1;
//
//
//                dividend.subtract(divider.clone().multiply(q, modulus), modulus)};
//
//                System.out.println(q);
//            }
//        }
////        System.out.println(invMod.getInverse(4));

    }

    @Test
    public void testXpowers() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 1);
        ArrayList<MutablePolynomial> exps = xPowers(poly, fastDivisionPreConditioning(poly, bigModulus), bigModulus);
        for (int i = 0; i < exps.size(); i++)
            assertTrue(exps.get(i).subtract(polyPowMod(createMonomial(1, bigModulus), i, poly, bigModulus, false)).isZero());
    }

    public static MutablePolynomial powerModWithXPowers(MutablePolynomial poly,
                                                        MutablePolynomial polyModulus,
                                                        InverseModMonomial invMod,
                                                        long modulus,
                                                        ArrayList<MutablePolynomial> xPowers) {
        poly = polyMod(poly, polyModulus, invMod, modulus, true);
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] == 0)
                continue;
            res.addMul(xPowers.get(i), poly.data[i], modulus);
        }
        return polyMod(res, polyModulus, invMod, modulus, false);
    }

    @Test
    public void testPowerMod0() throws Exception {
        int modulus = 19;
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        ArrayList<MutablePolynomial> exps = xPowers(polyModulus, invMod, modulus);
        MutablePolynomial poly = MutablePolynomial.create(3, 2, 12, 3, 4, 52, 12, 3, 1, 2, 3, 4);
        MutablePolynomial a = polyPowMod(poly, modulus, polyModulus, modulus, true);
        MutablePolynomial b = powerModWithXPowers(poly, polyModulus, invMod, modulus, exps);

        System.out.println(a);
        System.out.println(b);

        a = polyPowMod(a, modulus, polyModulus, modulus, true);
        b = powerModWithXPowers(b, polyModulus, invMod, modulus, exps);

        System.out.println(a);
        System.out.println(b);
    }

    public static ArrayList<MutablePolynomial> ExponentsNaive(MutablePolynomial polyModulus, long modulus, int n) {
        MutablePolynomial exponent = MutablePolynomial.create(0, 1);
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        for (int j = 0; j < n; j++) {
            exps.add(exponent);
            exponent = PolynomialArithmetics.polyPowMod(exponent, modulus, polyModulus, invMod, modulus, true);
        }
        return exps;
    }

    public static ArrayList<MutablePolynomial> ExponentsWithXPowers(MutablePolynomial polyModulus, long modulus, int n) {
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        return ExponentsWithXPowers(polyModulus, invMod, modulus, n, xPowers(polyModulus, invMod, modulus));
    }

    public static ArrayList<MutablePolynomial> ExponentsWithXPowers(MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, int n, ArrayList<MutablePolynomial> xPowers) {
        MutablePolynomial exponent = MutablePolynomial.create(0, 1);
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        for (int j = 0; j < n; j++) {
            exps.add(exponent);
            exponent = powerModWithXPowers(exponent, polyModulus, invMod, modulus, xPowers);
        }
        return exps;
    }

    public static ArrayList<MutablePolynomial> ExponentsWithComposition(MutablePolynomial polyModulus, long modulus, int n) {
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        MutablePolynomial base = polyMod(createMonomial(1, (int) modulus), polyModulus, invMod, modulus, false);
        MutablePolynomial exponent = base;
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        exps.add(MutablePolynomial.create(0, 1));
        for (int j = 0; ; j++) {
            exps.add(exponent);
            if (j == n - 1)
                break;
            exponent = horner(exponent, base, polyModulus, invMod, modulus);
        }
        return exps;
    }

    static MutablePolynomial horner(MutablePolynomial base, MutablePolynomial poly, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus) {
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = base.degree; i >= 0; --i)
//            res = LongArithmetics.add(LongArithmetics.multiply(res, point), data[i]);
            res = polyMod(res.multiply(poly, modulus).add(MutablePolynomial.create(base.data[i])), polyModulus, invMod, modulus, false);
        return res;
    }


    public static ArrayList<MutablePolynomial> hPowers(MutablePolynomial h, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, int nIterations) {
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        exps.add(MutablePolynomial.one());
        MutablePolynomial base = polyMod(h, polyModulus, invMod, modulus, true);
        exps.add(base);
        MutablePolynomial prev = base;
        for (int i = 0; i < nIterations; i++)
            exps.add(prev = polyMod(prev.clone().multiply(base, modulus), polyModulus, invMod, modulus, false));
        return exps;
    }

    static MutablePolynomial BrentKung(MutablePolynomial main, MutablePolynomial point, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus) {
        int t = (int) Math.sqrt(main.degree);
        ArrayList<MutablePolynomial> hPowers = hPowers(point, polyModulus, invMod, modulus, t);
        return BrentKung(main, point, polyModulus, invMod, modulus, t, hPowers);
    }


    static MutablePolynomial BrentKung(MutablePolynomial main, MutablePolynomial point, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, int t, ArrayList<MutablePolynomial> hPowers) {
        ArrayList<MutablePolynomial> gj = new ArrayList<>();
        for (int i = 0; i <= main.degree; ) {
            int to = i + t;
            if (to > (main.degree + 1))
                to = main.degree + 1;
            MutablePolynomial g = MutablePolynomial.create(Arrays.copyOfRange(main.data, i, to));
            gj.add(powerModWithXPowers(g, polyModulus, invMod, modulus, hPowers));
//            assert powerModWithXPowers(g, polyModulus, modulus, hPowers).subtract(horner(g, point, polyModulus, modulus), modulus).isZero();
            i = to;
        }
        MutablePolynomial pt = hPowers.get(t);
//        assert pt.clone().subtract(polyPowMod(point, t, polyModulus, modulus, true), modulus).isZero();
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = gj.size() - 1; i >= 0; --i)
            //            res = LongArithmetics.add(LongArithmetics.multiply(res, point), data[i]);
            res = polyMod(res.multiply(pt, modulus).add(gj.get(i)), polyModulus, invMod, modulus, false);
        return res;
    }

    public static ArrayList<MutablePolynomial> ExponentsWithBrentKung(MutablePolynomial polyModulus, long modulus, int n) {
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        MutablePolynomial base = polyMod(createMonomial(1, (int) modulus), polyModulus, invMod, modulus, false);
        MutablePolynomial exponent = base;
        ArrayList<MutablePolynomial> exps = new ArrayList<>();
        exps.add(MutablePolynomial.create(0, 1));
        int t = (int) Math.sqrt(polyModulus.degree);
        ArrayList<MutablePolynomial> hPowers = hPowers(base, polyModulus, invMod, modulus, t);
        for (int j = 0; ; j++) {
            exps.add(exponent);
            if (j == n - 1)
                break;
            exponent = BrentKung(exponent, base, polyModulus, invMod, modulus, t, hPowers);
        }
        return exps;
    }

    @Test
    public void test_h_powers() throws Exception {
        MutablePolynomial base = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 7, 8, 9);
        MutablePolynomial point = MutablePolynomial.create(1, 2, 3, 4);
        MutablePolynomial polyModulus = MutablePolynomial.create(4, 3, 2, 1);
        long modulus = 7;
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);

        ArrayList<MutablePolynomial> hPowers = hPowers(point, polyModulus, invMod, modulus, 4);
        System.out.println(hPowers.get(4).subtract(polyPowMod(point, 4, polyModulus, modulus, true)));
    }

    @Test
    public void test_brent_kung() throws Exception {
        MutablePolynomial base = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 7, 8, 9);
        MutablePolynomial point = MutablePolynomial.create(1, 2, 3, 4);
        MutablePolynomial polyModulus = MutablePolynomial.create(4, 3, 2, 1);
        long modulus = 7;

//        System.out.println(hPowers(polyModulus, point, modulus, 4));

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        System.out.println(horner(base, point, polyModulus, invMod, modulus));
        System.out.println("-----");
        System.out.println(BrentKung(base, point, polyModulus, invMod, modulus));
    }

    @Test
    public void test_brent_kung_2() throws Exception {
        MutablePolynomial base = bigPoly;
        MutablePolynomial point = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9).subtract(bigPoly, bigModulus);
        MutablePolynomial polyModulus = bigPoly;
        long modulus = bigModulus;

//        System.out.println(hPowers(polyModulus, point, modulus, 4));

        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        for (int i = 0; i < 1000; i++) {

            long start = System.nanoTime();
            MutablePolynomial a = horner(base, point, polyModulus, invMod, modulus);
            long horner = System.nanoTime() - start;
            start = System.nanoTime();
            MutablePolynomial b = BrentKung(base, point, polyModulus, invMod, modulus);
            long bk = System.nanoTime() - start;
            assertEquals(a, b);
            System.out.println("----");
            System.out.println(horner);
            System.out.println(bk);

        }
    }

    @Test
    public void exps_plain_vs_composition() throws Exception {
        ArrayList<MutablePolynomial> exponents = ExponentsNaive(bigPoly, bigModulus, 10);
        ArrayList<MutablePolynomial> exponentsx = ExponentsWithComposition(bigPoly, bigModulus, 10);
        ArrayList<MutablePolynomial> exponentsbk = ExponentsWithBrentKung(bigPoly, bigModulus, 10);
        for (int i = 0; i < exponents.size(); i++) {
            System.out.println(exponents.get(i));
            System.out.println(exponentsx.get(i));
            System.out.println(exponentsbk.get(i));
        }


    }

    @Test
    public void exps_plain_vs_xpowers_vs_horner_vs_brentkung() throws Exception {
        MutablePolynomial bigPoly = LargeDDFTest.bigPoly;
        MutablePolynomial factor = MutablePolynomial.create(4324, 2713, 2712, 4484, 1869, 1384, 452, 4955, 4523, 4247, 4261, 4549, 2584, 1);
        bigPoly = quotient(bigPoly, factor, bigModulus, true);
        for (int i = 0; i < 1000; i++) {

            System.out.println("------");

            long start = System.nanoTime();
            int nExponents = 120 - factor.degree;
            ArrayList<MutablePolynomial> exponents = ExponentsNaive(bigPoly, bigModulus, nExponents);
            System.out.println("naive:       " + (System.nanoTime() - start));


            start = System.nanoTime();
            InverseModMonomial invMod = fastDivisionPreConditioning(bigPoly, bigModulus);
            ArrayList<MutablePolynomial> xPowers = xPowers(bigPoly, invMod, bigModulus);
            ArrayList<MutablePolynomial> exponents_x = ExponentsWithXPowers(bigPoly, invMod, bigModulus, nExponents, xPowers);
            System.out.println("xPower:      " + (System.nanoTime() - start));


//            start = System.nanoTime();
//            ArrayList<MutablePolynomial> exponents_h = ExponentsWithComposition(bigPoly, bigModulus, nExponents);
//            System.out.println("composition: " + (System.nanoTime() - start));

            start = System.nanoTime();
            ArrayList<MutablePolynomial> exponents_bk = ExponentsWithBrentKung(bigPoly, bigModulus, nExponents);
            System.out.println("brend-kung:  " + (System.nanoTime() - start));

//            System.out.println(exponents.get(9));
//            System.out.println(exponents_x.get(9));
//            System.out.println(exponents_bk.get(9));
//            assertEquals(exponents.size(), exponents_x.size());
//            assertEquals(exponents.size(), exponents_h.size());
//            assertEquals(exponents.size(), exponents_bk.size());
        }
    }

    @Test
    public void test_tmp() throws Exception {
        ExponentsNaive(bigPoly, bigModulus, 20);
        System.out.println(1L * bigPoly.degree * bigPoly.degree * 20);

    }

    static long[][] qMatrixIterative(MutablePolynomial poly, int modulus) {
        int pDegree = poly.degree;
        long[][] matrix = new long[pDegree][pDegree];
        long[] prevRow = new long[pDegree], nextRow = new long[pDegree];
        prevRow[0] = 1;
        matrix[0][0] = 1;
        for (int i = 1, to = (int) ((pDegree - 1) * modulus); i <= to; ++i) {
            nextRow[0] = symMod(-prevRow[pDegree - 1] * poly.data[0], modulus);
            for (int j = 1; j < poly.degree; j++)
                nextRow[j] = symMod(prevRow[j - 1] - prevRow[pDegree - 1] * poly.data[j], modulus);
            if (i % modulus == 0)
                System.arraycopy(nextRow, 0, matrix[i / (int) modulus], 0, nextRow.length);
            long[] tmp = prevRow;
            prevRow = nextRow;
            nextRow = tmp;
        }
        return matrix;
    }


    static long[][] qMatrixHorner(MutablePolynomial poly, InverseModMonomial invMod, int modulus) {
        int pDegree = poly.degree;
        long[][] matrix = new long[pDegree][pDegree];
        matrix[0][0] = 1;
        MutablePolynomial xMod = polyMod(createMonomial(1, modulus), poly, modulus, false).symModulus(modulus);
        MutablePolynomial prev = xMod.clone();
        System.arraycopy(prev.data, 0, matrix[1], 0, pDegree);
        for (int i = 2; i < pDegree; ++i) {
            prev = horner(createMonomial(1, i), xMod, poly, invMod, modulus).symModulus(modulus);
            prev.ensureCapacity(pDegree);
            System.arraycopy(prev.data, 0, matrix[i], 0, pDegree);
        }
        return matrix;
    }


    static String toStringMatrix(long[][] matrix) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < matrix.length; i++) {
            sb.append(Arrays.toString(matrix[i])).append("\n");
        }
        return sb.toString();
    }

    @Test
    public void qma() throws Exception {
        MutablePolynomial xq = createMonomial(1, (int) bigModulus);
        xq = polyMod(xq, bigPoly, bigModulus, false);

        for (int i = 0; i < 10000; i++) {

            System.out.println("xxx");
            long start = System.nanoTime();
            InverseModMonomial invMod = fastDivisionPreConditioning(bigPoly, bigModulus);
            MutablePolynomial a = horner(xq, xq, bigPoly, invMod, bigModulus);
            System.out.println(System.nanoTime() - start);

            start = System.nanoTime();
            MutablePolynomial b = PolynomialArithmetics.polyPowMod(xq, bigModulus, bigPoly, bigModulus, true);
            System.out.println(System.nanoTime() - start);

            System.out.println(a.subtract(b));

        }
    }

    @Test
    public void name() throws Exception {
        long[][] expected = {
                {1, 0, 0, 0, 0, 0},
                {3, 5, -3, -3, -5, 5},
                {3, -5, -5, 1, -1, 0},
                {-2, 4, -1, 3, -4, -2},
                {-4, -3, -1, 0, 0, -3},
                {-3, -1, -4, -3, -1, -3}
        };
        MutablePolynomial poly = MutablePolynomial.create(1, -3, -1, -3, 1, -3, 1);
        int modulus = 11;
        long[][] qMatrix = qMatrixIterative(poly.clone(), modulus);
        System.out.println(toStringMatrix(qMatrix));
        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        System.out.println(toStringMatrix(qMatrixHorner(poly.clone(), invMod, modulus)));
        assertArrayEquals(expected, qMatrix);
    }
//
//    @Test
//    public void test1() throws Exception {
//        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 5);
//        MutablePolynomial x10 = createMonomial(1, 10);
//        MutablePolynomial x100 = createMonomial(1, 100);
//        MutablePolynomial x2 = createMonomial(1, 2);
//        MutablePolynomial x20 = createMonomial(1, 20);
//
//        assertEquals(
//                polyMod(x100, poly, bigModulus, true),
//                horner(x10, x10, poly, bigModulus)
//        );
//
//        assertEquals(
//                polyMod(x20, poly, bigModulus, true),
//                horner(x10, x2, poly, bigModulus)
//        );
//
////        if (true) return;
//
//        MutablePolynomial x10000 = createMonomial(1, 10000);
////        MutablePolynomial x100 = createMonomial(1, 100);
//        for (int i = 0; i < 10000; i++) {
//
//            System.out.println("====");
//            long start = System.nanoTime();
//            MutablePolynomial a = polyMod(x10000, bigPoly, bigModulus, true);
//            System.out.println(System.nanoTime() - start);
//            start = System.nanoTime();
//            MutablePolynomial b = horner(x100, x100, bigPoly, bigModulus);
//            System.out.println(System.nanoTime() - start);
//            Assert.assertEquals(a, b);
//        }
//    }
}
