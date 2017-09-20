package cc.redberry.rings;

/**
 * Defines {@link #parse(String)} method
 *
 * @since 1.0
 */
public interface ElementParser<E> {
    /**
     * Parse string into {@code E}
     */
    E parse(String string);
}
