package cc.redberry.rings.io;

/**
 * Defines {@link #parse(String)} method
 *
 * @since 1.0
 */
public interface IParser<Element> {
    /**
     * Parse string into {@code Element}
     */
    Element parse(String string);
}
