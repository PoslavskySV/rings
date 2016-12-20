package cc.r2.core.number;

public interface EuclideanRingElement<R extends EuclideanRingElement<R>>
        extends RingElement<R>, Comparable<R> {

    /**
     * Greatest common divisor of {@code this} and {@code oth}.
     *
     * @param oth other ring element
     * @return gcd of {@code this} and {@code oth}
     */
    R gcd(R oth);

    /**
     * Runs extended Euclidean algorithm to compute {@code [gcd(this, oth), x, y]} such that {@code x * this + y * oth = gcd(this, oth)}
     *
     * @param a a long
     * @param b a long
     * @return array of {@code [gcd(this, oth), x, y]} such that {@code x * this + y * oth = gcd(this, oth)}
     */
    R[] gcdExtended(R oth);

    /**
     * Least common multiple of{@code this} and {@code oth}.
     *
     * @param oth other ring element
     * @return gcd of {@code this} and {@code oth}
     */
    default R lcm(R oth) {
        R[] q = multiply(oth).divideAndRemainder(gcd(oth));
        assert q[1].isZero();
        return q[0];
    }

    /**
     * Returns an array of quotient and remainder two ring elements.
     *
     * @param oth other ring element
     */
    R[] divideAndRemainder(R oth);
}
