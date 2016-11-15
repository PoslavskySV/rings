package cc.r2.core.number;


public interface Ring<R extends RingElement> {
    R getOne();

    R getZero();

    R parse(Object o);

    Class<R> getElementType();
}
