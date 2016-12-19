package cc.r2.core.number;

abstract class ModBigIntegerRingAbstract<R extends RingElement> implements Ring<R> {
    final BigInteger mod;

    ModBigIntegerRingAbstract(BigInteger mod) {
        this.mod = mod;
    }

    public BigInteger getMod() {
        return mod;
    }

    @Override
    public R getOne() {
        ensureInitialized();
        return lazyHolder.ONE;
    }

    @Override
    public R getZero() {
        ensureInitialized();
        return lazyHolder.ZERO;
    }

    abstract R once_createOne();

    abstract R once_createZero();

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ModBigIntegerRingAbstract<?> that = (ModBigIntegerRingAbstract<?>) o;

        return mod.equals(that.mod);
    }

    @Override
    public int hashCode() {
        return mod.hashCode();
    }

    private volatile LazyHolder<R> lazyHolder = null;

    private void ensureInitialized() {
        if (lazyHolder == null)
            synchronized (this) {
                if (lazyHolder == null)
                    lazyHolder = new LazyHolder<>(once_createOne(), once_createZero());
            }
    }

    private static final class LazyHolder<R extends RingElement> {
        private final R ONE, ZERO;

        public LazyHolder(R ONE, R ZERO) {
            this.ONE = ONE;
            this.ZERO = ZERO;
        }
    }
}
