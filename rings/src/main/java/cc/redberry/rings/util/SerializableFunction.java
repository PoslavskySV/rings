package cc.redberry.rings.util;

/**
 *
 */
public interface SerializableFunction<T, R> extends java.io.Serializable {
    R apply(T t);

    default <V> SerializableFunction<T, V> andThen(SerializableFunction<? super R, ? extends V> after) {
        return (T t) -> after.apply(apply(t));
    }

    static <T> SerializableFunction<T, T> identity() {
        return t -> t;
    }
}
