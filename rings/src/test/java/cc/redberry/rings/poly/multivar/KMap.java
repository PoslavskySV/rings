package cc.redberry.rings.poly.multivar;

import com.koloboke.compile.KolobokeMap;

import java.util.Map;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
@KolobokeMap
public abstract class KMap implements Map<Long, MultivariatePolynomialZp64Test.PackedMonomial> {
    static KMap withExpectedSize(int expectedSize) {
        return new KolobokeKMap(expectedSize);
    }
}