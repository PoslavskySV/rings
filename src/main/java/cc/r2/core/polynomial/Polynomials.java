package cc.r2.core.polynomial;

import cc.r2.core.number.*;
import cc.r2.core.number.primes.SieveOfAtkin;
import gnu.trove.list.array.TIntArrayList;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.number.BigIntegerRing.IntegerRing;

/**
 * Created by poslavsky on 04/11/2016.
 */
public final class Polynomials {
    private Polynomials() {
    }

    @SuppressWarnings("unchecked")
    static <R extends RingElement<R>> R[] newArray(Ring<R> ring, int length) {
        return (R[]) Array.newInstance(ring.getElementType(), length);
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    R content(UnivariatePolynomial<R> poly) {
        if (poly.degree() == 0)
            return poly.coefficients[0];
        return cc.r2.core.number.Util.gcd(poly.coefficients);
    }

    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    primitivePart(UnivariatePolynomial<R> poly) {
        R content = content(poly);
        if (poly.lc().compareTo(poly.ring.getZero()) < 0)
            content = content.negate();
        return divide(poly, content);
    }

    private static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    divide(UnivariatePolynomial<R> poly, R content) {
        if (content.isOne())
            return poly;
        final R[] newData = newArray(poly.ring, poly.internalDegree);
        for (int i = 0; i < newData.length; i++) {
            R[] rs = poly.coefficients[i].divideAndRemainder(content);
            if (!rs[1].isZero())
                throw new IllegalArgumentException("" + rs[1]);
            newData[i] = rs[0];
        }
        return new UnivariatePolynomial<R>(poly.ring, newData);
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>[] pseudoDivideAndRemainder(
            final Ring<R> ring,
            UnivariatePolynomial<R> a,
            final UnivariatePolynomial<R> b) {
        return divideAndRemainder(ring, a.multiplyByFactor(b.lc().pow(a.degree() - b.degree() + 1)), b);
    }


    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>[] divideAndRemainder(
            final Ring<R> ring,
            final UnivariatePolynomial<R> a,
            final UnivariatePolynomial<R> b) {
        if (a.internalDegree < b.internalDegree)
            return null;

        MutableUnivariatePolynomial<R> r = new MutableUnivariatePolynomial<>(ring, a.degree(), a.coefficients.clone());

        final R[] q = newArray(ring, a.degree() - b.degree() + 1);
        for (int i = a.degree() - b.degree(); i >= 0; --i) {
            if (r.degree() == b.degree() + i) {
                final R[] qc = r.lc().divideAndRemainder(b.lc());
                if (!qc[1].isZero())
                    return null;
                q[i] = qc[0];
                r.subtract(b.coefficients, q[i], i);
            } else q[i] = ring.getZero();
        }
        return new UnivariatePolynomial[]{new UnivariatePolynomial<>(ring, q), r.toPoly()};
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> EuclideanPRS(Ring<R> ring,
                                               final UnivariatePolynomial<R> a,
                                               final UnivariatePolynomial<R> b) {
        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        UnivariatePolynomial<R> x = a, y = b, r;
        while (!y.isZero()) {
            UnivariatePolynomial<R>[] tmp = divideAndRemainder(ring, x, y);
            r = tmp[1];
            prs.add(r);
            x = y;
            y = r;
        }
        prs.add(x);
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> naivePRS(Ring<R> ring,
                                           final UnivariatePolynomial<R> a,
                                           final UnivariatePolynomial<R> b) {
        R aContent = content(a), bContent = content(b), contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        UnivariatePolynomial<R> x = aPP, y = bPP, r;
        while (!y.isZero()) {
            UnivariatePolynomial<R>[] tmp = pseudoDivideAndRemainder(ring, x, y);
            r = tmp[1];
            if (!r.isZero())
                r = primitivePart(r);
            prs.add(r);
            x = y;
            y = r;
        }
        if (x.degree() == 0)
            x = new UnivariatePolynomial<R>(ring, contentGCD);
        else
            x = primitivePart(x).multiplyByFactor(contentGCD);
        prs.add(x);
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> primitivePRS(final Ring<R> ring,
                                               final UnivariatePolynomial<R> a,
                                               final UnivariatePolynomial<R> b) {
        if (b.degree() > a.degree())
            return primitivePRS(ring, b, a);

        R aContent = content(a),
                bContent = content(b),
                contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        for (int i = 0; ; i++) {
            UnivariatePolynomial<R> ith = prs.get(i);
            UnivariatePolynomial<R> ithp1 = prs.get(i + 1);

            UnivariatePolynomial<R> q = divideAndRemainder(ring, ith.multiplyByFactor(ithp1.lc().pow(ith.degree() - ithp1.degree() + 1)), ithp1)[1];
            if (q.isZero())
                break;
            q = primitivePart(q);
            prs.add(q);
        }
        prs.set(prs.size() - 1, prs.get(prs.size() - 1).multiplyByFactor(contentGCD));
        return prs;
    }

    @SuppressWarnings({"unchecked", "ConstantConditions"})
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> subResultantPRS(final Ring<R> ring,
                                                  final UnivariatePolynomial<R> a,
                                                  final UnivariatePolynomial<R> b) {
        if (b.degree() > a.degree())
            return subResultantPRS(ring, b, a);

        R aContent = content(a),
                bContent = content(b),
                contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        ArrayList<R> beta = new ArrayList<>(), psi = new ArrayList<>();
        TIntArrayList deltas = new TIntArrayList();

        R cBeta, cPsi;
        for (int i = 0; ; i++) {
            UnivariatePolynomial<R> ith = prs.get(i);
            UnivariatePolynomial<R> ithp1 = prs.get(i + 1);
            int delta = ith.degree() - ithp1.degree();
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? ring.getOne() : ring.getOne().negate();
                cPsi = ring.getOne().negate();
            } else {
                cPsi = ith.lc().negate().pow(deltas.get(i - 1))
                ;//.multiply((R) BigInteger.valueOf(i + 1));
                if (deltas.get(i - 1) < 1)
                    cPsi = cPsi.multiply(psi.get(i - 1).pow(-deltas.get(i - 1) + 1));
                else {
                    R[] rs = cPsi.divideAndRemainder(psi.get(i - 1).pow(deltas.get(i - 1) - 1));
                    assert rs[1].isZero();
                    cPsi = rs[0];
                }
                cBeta = ith.lc().negate().multiply(cPsi.pow(delta));
            }

            UnivariatePolynomial<R> q = divideAndRemainder(ring, ith.multiplyByFactor(ithp1.lc().pow(ith.degree() - ithp1.degree() + 1)), ithp1)[1];
            if (q.isZero())
                break;

            q = divide(q, cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        prs.set(prs.size() - 1, prs.get(prs.size() - 1).multiplyByFactor(contentGCD));
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    EuclideanGCD(final Ring<R> ring,
                 final UnivariatePolynomial<R> a,
                 final UnivariatePolynomial<R> b) {
        if (b.internalDegree > a.internalDegree)
            return EuclideanGCD(ring, b, a);
        UnivariatePolynomial<R> x = a, y = b, r;
        while (!y.isZero()) {
            r = divideAndRemainder(ring, x, y)[1];
            x = y;
            y = r;
        }
        return x;
    }


    private static BigInteger bigIntSqRootCeil(BigInteger x)
            throws IllegalArgumentException {
        if (x.signum() < 0)
            throw new IllegalArgumentException("Negative argument.");

        if (x.isZero() || x.isOne())
            return x;
        BigInteger two = BigInteger.TWO;
        BigInteger y;
        for (y = x.divide(two);
             y.compareTo(x.divide(y)) > 0;
             y = ((x.divide(y)).add(y)).divide(two))
            ;
        if (x.compareTo(y.multiply(y)) == 0)
            return y;
        else
            return y.increment();
    }

    public static BigInteger norm(BigInteger[] factors) {
        BigInteger r = BigInteger.ZERO;
        for (BigInteger factor : factors)
            r = r.add(factor.multiply(factor));
        return bigIntSqRootCeil(r);
    }

    public static int log2(BigInteger f) {
        return f.bitLength();
    }

    public static BigInteger max(BigInteger a, BigInteger b) {
        return a.compareTo(b) < 0 ? b : a;
    }


    public static UnivariatePolynomial<ModPrimeBigInteger> convert(ModPrimeBigIntegerField primeF, UnivariatePolynomial<BigInteger> poly) {
        final ModPrimeBigInteger[] cfs = new ModPrimeBigInteger[poly.coefficients.length];
        for (int i = 0; i < cfs.length; i++)
            cfs[i] = new ModPrimeBigInteger(primeF, poly.coefficients[i]);
        return new UnivariatePolynomial<>(primeF, cfs);
    }

    public static UnivariatePolynomial<BigInteger> toSymMod(UnivariatePolynomial<ModPrimeBigInteger> poly) {
        final BigInteger[] cfs = new BigInteger[poly.coefficients.length];
        for (int i = 0; i < cfs.length; i++)
            cfs[i] = toSymMod(poly.coefficients[i].value(), poly.coefficients[i].getRing().getMod());
        return new UnivariatePolynomial<>(IntegerRing, cfs);
    }

    static BigInteger toSymMod(BigInteger b, BigInteger prime) {
//        BigInteger min = prime.divide(BigInteger.TWO);
//        BigInteger prev;
//        while (b.subtract(prime).compareTo(min) > 0){
//            prev = b;
//            b = b.subtract(prime);
//        }
//        return prev;
        BigInteger t;
        if(prime.mod(BigInteger.TWO).isZero())
            t = prime.divide(BigInteger.TWO);
        else
            t = prime.decrement().divide(BigInteger.TWO);
        if (b.compareTo(t) <= 0)
            return b;
        else return b.subtract(prime);
    }

    public static UnivariatePolynomial<BigInteger>
    modularGCD(final UnivariatePolynomial<BigInteger> f,
               final UnivariatePolynomial<BigInteger> g) {

        if (f.degree() < g.degree())
            return modularGCD(g, f);

        final BigInteger A = max(norm(f.coefficients), norm(g.coefficients));
        final int n = f.degree();
        final BigInteger bign = BigInteger.valueOf(n);
        final BigInteger b = f.lc().gcd(g.lc());
        final int k = 2 * log2(bign.increment().pow(n).multiply(b).multiply(A.pow(n * 2)));
        final BigInteger B = bigIntSqRootCeil(bign.increment()).multiply(BigInteger.TWO.pow(n)).multiply(A).multiply(b);
        final int l = log2(BigInteger.TWO.multiply(B).increment());


        ArrayList<UnivariatePolynomial<ModPrimeBigInteger>> mods = new ArrayList<>();
        ArrayList<BigInteger> primes = new ArrayList<>();

        int TTT = (int) (2 * k * Math.log(k));
        SieveOfAtkin sieveOfAtkin = SieveOfAtkin.createSieve(TTT);
        int minDegree = -1;
        for (int i = 0; i < TTT; i++) {
            if (!sieveOfAtkin.isPrime(i))
                continue;
            BigInteger prime = BigInteger.valueOf(i);
//            System.out.println("->" + prime);
            if (b.divideAndRemainder(prime)[1].isZero())
                continue;


            ModPrimeBigIntegerField field = new ModPrimeBigIntegerField(prime);
            UnivariatePolynomial<ModPrimeBigInteger> fMod = convert(field, f);
            UnivariatePolynomial<ModPrimeBigInteger> gMod = convert(field, g);
            UnivariatePolynomial<ModPrimeBigInteger> modGcd = EuclideanGCD(field, fMod, gMod);

            if (modGcd.degree() == 0)
                return f.getOne();


            if (minDegree == -1)
                minDegree = modGcd.degree();
            else if (modGcd.degree() < minDegree) {
                minDegree = modGcd.degree();
                mods.clear();
                primes.clear();
            } else if (modGcd.degree() > minDegree)
                continue;



            UnivariatePolynomial<ModPrimeBigInteger> tmoMOd = modGcd.multiplyByFactor(new ModPrimeBigInteger(field, b));
            mods.add(tmoMOd);
            primes.add(prime);

            System.out.println(prime);
            System.out.println(tmoMOd);
            System.out.println("====");

            if (mods.size() > l) {
                //CRT
                BigInteger[] newCoeffs = new BigInteger[minDegree + 1];
                BigInteger[] coprimes = primes.toArray(new BigInteger[primes.size()]);
                BigInteger hui = BigInteger.ONE;
                for (int j = 0; j < coprimes.length; j++)
                    hui = hui.multiply(coprimes[j]);

                System.out.println("HUI: " + hui);
                BigInteger[] tmpCfx = new BigInteger[primes.size()];
                for (int j = 0; j < minDegree + 1; j++) {
                    for (int m = 0; m < tmpCfx.length; m++)
                        tmpCfx[m] = mods.get(m).coefficients[j].value();

                    System.out.println(Arrays.toString(coprimes));
                    System.out.println(Arrays.toString(tmpCfx));
                    BigInteger tmp = ChineseRemainderAlgorithm.CRT(coprimes, tmpCfx).mod(hui);
                    System.out.println(tmp);
                    newCoeffs[j] = toSymMod(tmp,hui);
                    System.out.println(newCoeffs[j]);
                }


                UnivariatePolynomial<BigInteger> poly = new UnivariatePolynomial<>(IntegerRing, newCoeffs);

                System.out.println("----");
                System.out.println(hui);
                System.out.println("----");
                System.out.println(poly);
                System.out.println("----");
                System.out.println(primitivePart(poly));
                if (divideAndRemainder(IntegerRing, f, poly)[1].isZero() && divideAndRemainder(IntegerRing, g, poly)[1].isZero())
                    return poly;
            }
        }
//        BigInteger k = BigInteger.TWO.multiply(n.increment().pow(n).multiply(b).multiply());
        return null;
    }
}
