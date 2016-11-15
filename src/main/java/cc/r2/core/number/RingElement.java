package cc.r2.core.number;


public interface RingElement<R extends RingElement<R>> {
    R add(R a);

    R subtract(R a);

    R multiply(R a);

    @SuppressWarnings("unchecked")
    default R negate() {
        return getZero().subtract((R) this);
    }

    R getZero();

    R getOne();

    default R increment() {
        return add(getOne());
    }

    default R decrement() {
        return subtract(getOne());
    }

    @SuppressWarnings("unchecked")
    default R pow(int p) {
        if (p < 0)
            throw new IllegalArgumentException();
        if (p == 0)
            return getOne();
        if (p == 1)
            return (R) this;

        R result = getOne();
        R k2p = (R) this;
        while (p != 0) {
            if ((p & 0x1) != 0)
                result = result.multiply(k2p);
            k2p = k2p.multiply(k2p);
            p >>= 1;
        }
        return result;
    }

    Ring<R> getRing();

    default boolean isZero() {
        return this.equals(getZero());
    }

    default boolean isOne() {
        return this.equals(getOne());
    }
}
