package cc.r2.core.polynomial;

import cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.InverseModMonomial;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.MutableLongPoly.createMonomial;
import static cc.r2.core.polynomial.SmallPolynomialArithmetics.*;
import static cc.r2.core.polynomial.SmallPolynomials.DistinctDegreeFactorization;
import static cc.r2.core.polynomial.SmallPolynomials.PolynomialGCD;
import static cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.quotient;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 16/01/2017.
 */
public class LargeDDFtest2 {
    static final MutableLongPoly bigPoly = MutableLongPoly.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1);
    static final int bigModulus = 5659;


    /** returns x^{i*modulus} mod polyModulus for i in [0...degree] */
    public static ArrayList<MutableLongPoly> xPowers(MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus) {
        return hPowers(createMonomial(1, (int) modulus), polyModulus, invMod, modulus, polyModulus.degree);
    }

    /** returns h^{i*modulus} mod polyModulus for i in [0...nIterations] */
    public static ArrayList<MutableLongPoly> hPowers(MutableLongPoly h, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, int nIterations) {
        ArrayList<MutableLongPoly> exponents = new ArrayList<>();
        hPowers(h, polyModulus, invMod, modulus, nIterations, exponents);
        return exponents;
    }

    /** writes h^{i*modulus} mod polyModulus for i in [0...nIterations] to exponents */
    public static void hPowers(MutableLongPoly h, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, int nIterations, ArrayList<MutableLongPoly> exponents) {
        exponents.add(MutableLongPoly.one());
        MutableLongPoly base = polyMod(h, polyModulus, invMod, modulus, true);
        exponents.add(base);
        MutableLongPoly prev = base;
        for (int i = 0; i < nIterations; i++)
            exponents.add(prev = polyMod(prev.clone().multiply(base, modulus), polyModulus, invMod, modulus, false));
    }


    @Test
    public void testXPowers() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 2, 3, 4, 5, 6, 1);
        ArrayList<MutableLongPoly> exps = xPowers(poly, fastDivisionPreConditioning(poly, bigModulus), bigModulus);
        for (int i = 0; i < exps.size(); i++)
            assertTrue(exps.get(i).subtract(polyPowMod(createMonomial(1, bigModulus), i, poly, bigModulus, false)).isZero());
    }


    /** returns {@code poly^modulus mod polyModulus} **/
    public static MutableLongPoly raiseWithXPowers(MutableLongPoly poly,
                                                   MutableLongPoly polyModulus,
                                                   InverseModMonomial invMod,
                                                   long modulus,
                                                   ArrayList<MutableLongPoly> xPowers) {
        poly = polyMod(poly, polyModulus, invMod, modulus, true);
        MutableLongPoly res = MutableLongPoly.zero();
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] == 0)
                continue;
            res.addMul(xPowers.get(i), poly.data[i], modulus);
        }
        return polyMod(res, polyModulus, invMod, modulus, false);
    }

    @Test
    public void testRaiseWithXPowers() throws Exception {
        MutableLongPoly poly = bigPoly;
        long modulus = bigModulus;
        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        ArrayList<MutableLongPoly> xPowers = xPowers(poly, invMod, modulus);

        MutableLongPoly raised = poly;
        for (int i = 0; i < 15; i++) {
            MutableLongPoly prev = raised;
            MutableLongPoly next = raiseWithXPowers(raised, poly, invMod, modulus, xPowers);
            assertEquals(next, polyPowMod(prev, modulus, poly, invMod, modulus, true));
            raised = next;
        }


    }

    /** returns main(point) mod polyModulus */
    static MutableLongPoly BrentKung(MutableLongPoly main, MutableLongPoly point, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus) {
        int t = (int) Math.sqrt(main.degree);
        ArrayList<MutableLongPoly> hPowers = hPowers(point, polyModulus, invMod, modulus, t);
        return BrentKung(main, polyModulus, invMod, modulus, t, hPowers);
    }

    /** returns main(point) mod polyModulus */
    static MutableLongPoly BrentKung(
            MutableLongPoly main,
            MutableLongPoly polyModulus,
            InverseModMonomial invMod,
            long modulus,
            int tBrentKung,
            ArrayList<MutableLongPoly> pointPowers) {
        ArrayList<MutableLongPoly> gj = new ArrayList<>();
        for (int i = 0; i <= main.degree; ) {
            int to = i + tBrentKung;
            if (to > (main.degree + 1))
                to = main.degree + 1;
            MutableLongPoly g = MutableLongPoly.create(Arrays.copyOfRange(main.data, i, to));
            gj.add(raiseWithXPowers(g, polyModulus, invMod, modulus, pointPowers));
            i = to;
        }
        MutableLongPoly pt = pointPowers.get(tBrentKung);
        MutableLongPoly res = MutableLongPoly.zero();
        for (int i = gj.size() - 1; i >= 0; --i)
            res = polyMod(res.multiply(pt, modulus).add(gj.get(i)), polyModulus, invMod, modulus, false);
        return res;
    }

    static final class GiantSteps {
        final ArrayList<MutableLongPoly> giantSteps = new ArrayList<>();
        final MutableLongPoly poly;
        final InverseModMonomial invMod;
        final long modulus;
        final MutableLongPoly basePower;

        final int tBrentKung;
        final ArrayList<MutableLongPoly> hPowers = new ArrayList<>();

        public GiantSteps(MutableLongPoly poly, InverseModMonomial invMod, long modulus, MutableLongPoly basePower) {
            this.poly = poly;
            this.invMod = invMod;
            this.modulus = modulus;
            this.basePower = basePower;
            this.tBrentKung = (int) Math.sqrt(poly.degree);

            giantSteps.add(MutableLongPoly.createMonomial(1, 1)); // <- add x
            giantSteps.add(basePower); // <- add x^p mod poly
        }

        MutableLongPoly get(int j) {
            if (j < giantSteps.size())
                return giantSteps.get(j);
            if (hPowers.isEmpty())
                hPowers(basePower, poly, invMod, modulus, tBrentKung, hPowers);

            MutableLongPoly xPowerBig = giantSteps.get(giantSteps.size() - 1);
            for (int i = giantSteps.size(); i <= j; ++i)
                giantSteps.add(xPowerBig = BrentKung(xPowerBig, poly, invMod, modulus, tBrentKung, hPowers));
            return xPowerBig;
        }
    }

    static final class DDFSteps2 {
        final int B, l, m;
        final ArrayList<MutableLongPoly> babySteps;
        final GiantSteps giantSteps;
        final InverseModMonomial invMod;

        public DDFSteps2(int b, int l, int m, ArrayList<MutableLongPoly> babySteps, GiantSteps giantSteps, InverseModMonomial invMod) {
            B = b;
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    static DDFSteps2 ShoupDDF2(MutableLongPoly poly, long modulus) {
        int n = poly.degree;
        int B = (int) Math.floor(n / 2.);
        int l = (int) Math.floor(Math.sqrt(B));
        int m = (int) Math.ceil(1.0 * B / l);

        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        ArrayList<MutableLongPoly> xPowers = xPowers(poly, invMod, modulus);

        //baby steps
        ArrayList<MutableLongPoly> babySteps = new ArrayList<>();
        babySteps.add(MutableLongPoly.createMonomial(1, 1)); // <- add x
        MutableLongPoly xPower = xPowers.get(1); // x^p mod poly
        babySteps.add(xPower); // <- add x^p mod poly
        for (int i = 0; i <= l - 2; ++i)
            babySteps.add(xPower = raiseWithXPowers(xPower, poly, invMod, modulus, xPowers));

        // <- xPower = x^(p^l) mod poly
        return new DDFSteps2(B, l, m, babySteps, new GiantSteps(poly, invMod, modulus, xPower), invMod);
    }

    static final class IBases {
        final DDFSteps2 steps;
        final long modulus;
        final MutableLongPoly poly;

        public IBases(DDFSteps2 steps, long modulus, MutableLongPoly poly) {
            this.steps = steps;
            this.modulus = modulus;
            this.poly = poly;
        }

        final ArrayList<MutableLongPoly> iBases = new ArrayList<>();

        public MutableLongPoly get(int k) {
            if (k < iBases.size())
                return iBases.get(k);

            for (int j = iBases.size(); j <= k; ++j) {
                MutableLongPoly iBase = MutableLongPoly.one();
                for (int i = 0; i <= steps.l - 1; ++i) {
                    MutableLongPoly tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus);
                    iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, modulus, false);
                }
                iBases.add(iBase);
            }
            return iBases.get(k);
        }
    }

    static ArrayList<MutableLongPoly> ShoupDDFFactorization2(MutableLongPoly poly, long modulus) {
        return ShoupDDFFactorization2(poly, modulus, ShoupDDF2(poly, modulus));
    }


    private static ArrayList<MutableLongPoly> ShoupDDFFactorization2(MutableLongPoly poly, long modulus, DDFSteps2 steps) {
        IBases iBases = new IBases(steps, modulus, poly);

        ArrayList<MutableLongPoly> fList = new ArrayList<>();
        for (int i = 0; i <= poly.degree; ++i)
            fList.add(MutableLongPoly.one());

        MutableLongPoly current = poly.clone();
        for (int j = 1; j <= steps.m; ++j) {
            MutableLongPoly gcd = PolynomialGCD(current, iBases.get(j), modulus);
            if (gcd.isConstant())
                continue;
            current = quotient(current, gcd, modulus, false);
            for (int i = steps.l - 1; i >= 0; --i) {
                MutableLongPoly tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus), modulus);
                fList.set(steps.l * j - i, tmp);
                gcd = quotient(gcd, tmp, modulus, false);
            }
        }
        if (!current.isOne())
            fList.set(current.degree - 1, current);

        for (int i = fList.size() - 1; i >= 0; --i)
            if (fList.get(i).isOne())
                fList.remove(i);

        return fList;
    }

    @Test
    public void testGiantSteps() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 2, 3, 4, 5, 5, 6, 5, 4, 3, 2, 1);
        long modulus = 3;
        int l = 3;
        MutableLongPoly x = MutableLongPoly.createMonomial(1, 1);
        MutableLongPoly basePower = polyPowMod(x, (int) LongArithmetics.pow(modulus, l), poly, modulus, true);
        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        GiantSteps giantSteps = new GiantSteps(poly, invMod, modulus, basePower);

        for (int j = 0; j < 7; j++) {
            MutableLongPoly expected = polyPowMod(x, LongArithmetics.toInt(LongArithmetics.pow(modulus, l * j)), poly, modulus, true);
            assertEquals(expected, giantSteps.get(j));
        }
    }

    @Test
    public void testDDF3_2() throws Exception {
        Well1024a rnd = new Well1024a();
        long modulus = 17;
        for (int i = 0; i < 1000; i++) {
            System.out.println("====");
            MutableLongPoly poly = RandomPolynomials.randomMonicPoly(250, modulus, rnd);
            long start = System.nanoTime();
            ArrayList<MutableLongPoly> f1 = ShoupDDFFactorization(poly, modulus);
            System.out.println(System.nanoTime() - start);

            start = System.nanoTime();
            ArrayList<MutableLongPoly> f2 = ShoupDDFFactorization2(poly, modulus);
            System.out.println(System.nanoTime() - start);

            start = System.nanoTime();
            DistinctDegreeFactorization(poly, modulus);
            System.out.println(System.nanoTime() - start);

            assertEquals(f1, f2);
        }
//        long modulus = bigModulus;
//        MutableLongPoly poly = bigPoly;//.clone().multiply(bigPoly, modulus).square(modulus).add(bigPoly.cut(10), modulus);
//        poly.modulus(modulus).monic(modulus);
//        System.out.println(poly);
//        System.out.println(poly.degree);
//
//        ArrayList<MutableLongPoly> f = ShoupDDFFactorization2(poly, modulus);
//        System.out.println(f.size());
//        MutableLongPoly e = MutableLongPoly.one();
//        for (MutableLongPoly factor : f) {
//            System.out.println(factor);
//            e.multiply(factor, modulus);
//        }
//        System.out.println(e.equals(poly));
    }

    static DDFSteps ShoupDDF(MutableLongPoly poly, long modulus) {
        int n = poly.degree;
        int B = (int) Math.floor(n / 2.);
        int l = (int) Math.floor(Math.sqrt(B));
        int m = (int) Math.ceil(1.0 * B / l);

        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);
        ArrayList<MutableLongPoly> xPowers = xPowers(poly, invMod, modulus);

        //baby steps
        ArrayList<MutableLongPoly> babySteps = new ArrayList<>();
        babySteps.add(MutableLongPoly.createMonomial(1, 1)); // <- add x
        MutableLongPoly xPower = xPowers.get(1); // x^p mod poly
        babySteps.add(xPower); // <- add x^p mod poly
        for (int i = 0; i <= l - 2; ++i)
            babySteps.add(xPower = raiseWithXPowers(xPower, poly, invMod, modulus, xPowers));

        // <- xPower = x^(p^l) mod poly

        //giant steps
        ArrayList<MutableLongPoly> giantSteps = new ArrayList<>();
        giantSteps.add(MutableLongPoly.createMonomial(1, 1)); // <- add x
        giantSteps.add(xPower);
        MutableLongPoly xPowerBig = xPower;
        int tBrentKung = (int) Math.sqrt(poly.degree);
        ArrayList<MutableLongPoly> hPowers = hPowers(xPowerBig, poly, invMod, modulus, tBrentKung);
        for (int i = 0; i < m - 1; ++i)
            giantSteps.add(xPowerBig = BrentKung(xPowerBig, poly, invMod, modulus, tBrentKung, hPowers));

        return new DDFSteps(B, l, m, babySteps, giantSteps, invMod);
    }

    static ArrayList<MutableLongPoly> ShoupDDFFactorization(MutableLongPoly poly, long modulus) {
        return ShoupDDFFactorization(poly, modulus, ShoupDDF(poly, modulus));
    }

    private static ArrayList<MutableLongPoly> ShoupDDFFactorization(MutableLongPoly poly, long modulus, DDFSteps steps) {
        ArrayList<MutableLongPoly> iBases = new ArrayList<>();
        for (int j = 0; j <= steps.m; ++j) {
            MutableLongPoly iBase = MutableLongPoly.one();
            for (int i = 0; i <= steps.l - 1; ++i) {
                MutableLongPoly tmp = steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus);
                iBase = polyMultiplyMod(iBase, tmp, poly, steps.invMod, modulus, false);
            }
            iBases.add(iBase);
        }

        ArrayList<MutableLongPoly> fList = new ArrayList<>();
        for (int i = 0; i <= poly.degree; ++i)
            fList.add(MutableLongPoly.one());

        MutableLongPoly current = poly.clone();
        for (int j = 1; j <= steps.m; ++j) {
            MutableLongPoly gcd = PolynomialGCD(current, iBases.get(j), modulus);
            if (gcd.isConstant())
                continue;
            current = quotient(current, gcd, modulus, false);
            for (int i = steps.l - 1; i >= 0; --i) {
                MutableLongPoly tmp = PolynomialGCD(gcd, steps.giantSteps.get(j).clone().subtract(steps.babySteps.get(i), modulus), modulus);
                fList.set(steps.l * j - i, tmp);
                gcd = quotient(gcd, tmp, modulus, false);
            }
        }
        if (!current.isOne())
            fList.set(current.degree - 1, current);

        for (int i = fList.size() - 1; i >= 0; --i)
            if (fList.get(i).isOne())
                fList.remove(i);


        return fList;
    }


    static final class DDFSteps {
        final int B, l, m;
        final ArrayList<MutableLongPoly> babySteps;
        final ArrayList<MutableLongPoly> giantSteps;
        final InverseModMonomial invMod;

        public DDFSteps(int b, int l, int m, ArrayList<MutableLongPoly> babySteps, ArrayList<MutableLongPoly> giantSteps, InverseModMonomial invMod) {
            B = b;
            this.l = l;
            this.m = m;
            this.babySteps = babySteps;
            this.giantSteps = giantSteps;
            this.invMod = invMod;
        }
    }

    @Test
    public void testDDF1() throws Exception {
        long modulus = 3;
        MutableLongPoly poly = RandomPolynomials.randomMonicPoly(50, modulus, new Well1024a());
        DDFSteps ddfSteps = ShoupDDF(poly, modulus);

        for (int i = 0; i <= ddfSteps.l; i++)
            assertEquals(polyMod(
                    MutableLongPoly.createMonomial(1, (int) LongArithmetics.pow(modulus, i)), poly, modulus, false), ddfSteps.babySteps.get(i));

        for (int i = 0; i <= ddfSteps.m; i++) {
            MutableLongPoly expected = polyPowMod(MutableLongPoly.createMonomial(1, 1), LongArithmetics.pow(modulus, ddfSteps.l * i), poly, modulus, false);
            assertEquals(expected, ddfSteps.giantSteps.get(i));
        }
    }

    @Test
    public void testDDF1a() throws Exception {
        long modulus = 2;
        MutableLongPoly poly = bigPoly.clone().cut(50);
        poly.modulus(modulus).monic(modulus);
        DDFSteps ddfSteps = ShoupDDF(poly, modulus);

        for (int i = 0; i <= ddfSteps.l; i++) {
            MutableLongPoly expected = polyPowMod(MutableLongPoly.createMonomial(1, 1), LongArithmetics.pow(modulus, i), poly, modulus, false);
            assertEquals(expected, ddfSteps.babySteps.get(i));
        }

        for (int i = 0; i <= ddfSteps.m; i++) {
            MutableLongPoly expected = polyPowMod(MutableLongPoly.createMonomial(1, 1), LongArithmetics.pow(modulus, ddfSteps.l * i), poly, modulus, false);
            assertEquals(expected, ddfSteps.giantSteps.get(i));
        }
    }

    @Test
    public void testDDF2() throws Exception {
        for (int i = 0; i < 100; i++) {
            long start = System.nanoTime();
            ShoupDDFFactorization(bigPoly, bigModulus);
            System.out.println(System.nanoTime() - start);
        }
    }

    @Test
    public void testDDF3() throws Exception {
        long modulus = bigModulus;
//        System.out.println(bigPoly);
//        MutableLongPoly poly = RandomPolynomials.randomMonicPoly(12, bigModulus, new Well1024a()).multiply(RandomPolynomials.randomMonicPoly(250,modulus,new Well1024a()));
        MutableLongPoly poly = bigPoly;//.clone().multiply(bigPoly, modulus).square(modulus).add(bigPoly.cut(10), modulus);
//        poly = polyPowMod(poly, 10, modulus, false).add(poly.clone().cut(100), modulus).add(MutableLongPoly.one(), modulus);
//        poly = poly.cut(3000);
        poly.modulus(modulus).monic(modulus);
        System.out.println(poly);
        System.out.println(poly.degree);

        DDFSteps ddfSteps = ShoupDDF(poly, modulus);
        ArrayList<MutableLongPoly> f = ShoupDDFFactorization(poly, modulus, ddfSteps);
        System.out.println(f.size());
//        System.out.println(f);
        MutableLongPoly e = MutableLongPoly.one();
        for (MutableLongPoly factor : f) {
            System.out.println(factor);
            e.multiply(factor, modulus);
        }
        System.out.println(e.equals(poly));
    }
}
