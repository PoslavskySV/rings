package cc.redberry.rings.poly;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.univar.*;

import java.util.Arrays;
import java.util.Iterator;

/**
 * Galois field {@code GF(p, q)}. Galois field is represented as a quotient ring {@code F[x]/<m(x)>} with given
 * irreducible polynomial {@code m(x)} (minimal polynomial); cardinality of then field is than {@code p^q} where {@code
 * p} is the cardinality of {@code F} and {@code q} is the degree of minimal polynomial. Since Galois field is in fact a
 * simple field extension it inherits all corresponding methods.
 *
 * @param <E> type of polynomials that represent elements of this Galois field
 * @see AlgebraicNumberField
 * @see cc.redberry.rings.Rings#GF(IUnivariatePolynomial)
 * @see cc.redberry.rings.Rings#GF(long, int)
 * @since 1.0
 */
public final class FiniteField<E extends IUnivariatePolynomial<E>>
        extends SimpleFieldExtension<E> {
    private static final long serialVersionUID = 1L;
    /** GF(3^3) */
    public static final FiniteField<UnivariatePolynomialZp64> GF27 = new FiniteField<>(UnivariatePolynomialZ64.create(-1, -1, 0, 1).modulus(3));
    /** GF(17^5) */
    public static final FiniteField<UnivariatePolynomialZp64> GF17p5 = new FiniteField<>(UnivariatePolynomialZ64.create(11, 11, 0, 3, 9, 9).modulus(17).monic());

    /**
     * Constructs finite field from the specified irreducible polynomial.
     *
     * <p><b>NOTE:</b> irreducibility test for the minimal polynomial is not performed here, use {@link
     * IrreduciblePolynomials#irreducibleQ(IUnivariatePolynomial)} to test irreducibility.
     *
     * @param minimalPoly the minimal polynomial
     */
    public FiniteField(E minimalPoly) {
        super(minimalPoly);
        if (!minimalPoly.isOverFiniteField())
            throw new IllegalArgumentException("Irreducible poly must be over finite field.");
    }

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public boolean isUnit(E element) {
        return !element.isZero();
    }

    @Override
    public E gcd(E a, E b) {
        return a;
    }

    @Override
    public E[] divideAndRemainder(E a, E b) {
        return a.createArray(multiply(a, reciprocal(b)), getZero());
    }

    @Override
    public E remainder(E dividend, E divider) {
        return getZero();
    }

    /**
     * Returns iterator over all field elements
     *
     * @return iterator over all field elements
     */
    @SuppressWarnings("unchecked")
    @Override
    public Iterator<E> iterator() {
        if (!isFinite())
            throw new RuntimeException("Ring of infinite cardinality.");
        if (minimalPoly instanceof UnivariatePolynomial)
            return (Iterator<E>) new It(((UnivariatePolynomial) minimalPoly).ring, minimalPoly.degree());
        else if (minimalPoly instanceof UnivariatePolynomialZp64)
            return (Iterator<E>) new lIt(((UnivariatePolynomialZp64) minimalPoly).ring, minimalPoly.degree());
        throw new RuntimeException();
    }

    private static final class It<E> implements Iterator<UnivariatePolynomial<E>> {
        final Ring<E> ring;
        final E[] data;
        final Iterator<E>[] iterators;

        @SuppressWarnings("unchecked")
        It(Ring<E> ring, int degree) {
            this.ring = ring;
            this.data = ring.createArray(degree);
            this.iterators = new Iterator[degree];
            for (int i = 0; i < iterators.length; i++)
                iterators[i] = ring.iterator();
            for (int i = 0; i < data.length; i++)
                data[i] = iterators[i].next();
        }

        @Override
        public boolean hasNext() {
            return Arrays.stream(iterators).anyMatch(Iterator::hasNext);
        }

        private boolean first = true;

        @Override
        public UnivariatePolynomial<E> next() {
            if (first) {
                first = false;
                return UnivariatePolynomial.create(ring, data.clone());
            }
            int i = 0;
            if (!iterators[i].hasNext())
                while (i < iterators.length && !iterators[i].hasNext()) {
                    iterators[i] = ring.iterator();
                    data[i] = iterators[i].next();
                    ++i;
                }

            if (i >= iterators.length)
                return null;

            data[i] = iterators[i].next();
            return UnivariatePolynomial.createUnsafe(ring, data.clone());
        }
    }

    private static final class lIt implements Iterator<UnivariatePolynomialZp64> {
        final IntegersZp64 ring;
        final long[] data;

        @SuppressWarnings("unchecked")
        lIt(IntegersZp64 ring, int degree) {
            this.ring = ring;
            this.data = new long[degree];
        }

        @Override
        public boolean hasNext() {
            return Arrays.stream(data).anyMatch(l -> l < (ring.modulus - 1));
        }

        private boolean first = true;

        @Override
        public UnivariatePolynomialZp64 next() {
            if (first) {
                first = false;
                return UnivariatePolynomialZp64.createUnsafe(ring, data.clone());
            }
            int i = 0;
            if (data[i] >= ring.modulus - 1)
                while (i < data.length && data[i] >= ring.modulus - 1) {
                    data[i] = 0;
                    ++i;
                }

            if (i >= data.length)
                return null;

            ++data[i];
            return UnivariatePolynomialZp64.createUnsafe(ring, data.clone());
        }
    }
}
