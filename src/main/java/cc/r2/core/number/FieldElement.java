package cc.r2.core.number;


public interface FieldElement<R extends FieldElement<R>>
        extends EuclideanRingElement<R> {
    R divide(R a);

    @SuppressWarnings("unchecked")
    default R reciprocal() {
        return getOne().divide((R) this);
    }
}
