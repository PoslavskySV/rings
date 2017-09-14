package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface ElementParser<E> {
    E parse(String string);
}
