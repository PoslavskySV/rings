package cc.r2.core.combinatorics;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class InconsistentGeneratorsException extends RuntimeException {
    public InconsistentGeneratorsException(String message) {
        super(message);
    }

    public InconsistentGeneratorsException() {
    }
}
