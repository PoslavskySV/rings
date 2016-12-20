package cc.r2.core.number;

import org.apache.commons.math3.util.ArithmeticUtils;

public final class Util {
    private Util() {
    }


    public static <R extends EuclideanRingElement<R>> R gcd(final R... vals) {
        if (vals.length < 2)
            throw new IllegalArgumentException();
        R gcd = vals[0].gcd(vals[1]);
        for (int i = 2; i < vals.length; i++)
            gcd = vals[i].gcd(gcd);
        return gcd;
    }

    public static <R extends EuclideanRingElement<R>> R lcm(final R... vals) {
        if (vals.length < 2)
            throw new IllegalArgumentException();
        R lcm = vals[0].lcm(vals[1]);
        for (int i = 2; i < vals.length; i++)
            lcm = vals[i].lcm(lcm);
        return lcm;
    }
}
