package cc.r2.core.number;

/**
 * Created by poslavsky on 06/12/2016.
 */
public final class ChineseRemainderAlgorithm {
    public static BigInteger CRT(final BigInteger[] coprimes, final BigInteger[] remainders) {
        if (coprimes.length != remainders.length)
            throw new IllegalArgumentException();
        BigInteger m = coprimes[0];
        for (int i = 1; i < coprimes.length; i++) {
            if (coprimes[i].signum() <= 0)
                throw new RuntimeException("Negative CRT input: " + coprimes[i]);
            m = coprimes[i].multiply(m);
        }

        BigInteger result = BigInteger.ZERO;
        for (int i = 0; i < coprimes.length; i++) {
            BigInteger mi = m.divide(coprimes[i]);
            BigInteger eea = bezout0(mi, coprimes[i]);
            result = result.add(mi.multiply(eea.multiply(remainders[i]).mod(coprimes[i])));
        }
        return result;
    }

    private static BigInteger bezout0(BigInteger a, BigInteger b) {
        BigInteger s = BigInteger.ZERO, old_s = BigInteger.ONE;
        BigInteger r = b, old_r = a;

        BigInteger q;
        BigInteger tmp;
        while (!r.isZero()) {
            q = old_r.divide(r);

            tmp = old_r;
            old_r = r;
            r = tmp.subtract(q.multiply(r));

            tmp = old_s;
            old_s = s;
            s = tmp.subtract(q.multiply(s));
        }
        assert old_r.isOne();
        return old_s;
    }
}
