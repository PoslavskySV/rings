package cc.redberry.rings.io;

/**
 * Elements that could be stringified with the help of IStringifier
 */
public interface Stringifiable<E> {
    /** convert this to string with the use of stringifier */
    String toString(IStringifier<E> stringifier);
}
