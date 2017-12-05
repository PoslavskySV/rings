package cc.redberry.rings;

import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Factor decomposition of element. Unit coefficient of decomposition is stored in {@link #unit}, factors returned by
 * {@link #get(int)} are non-units. This class is mutable.
 *
 * <p><i>Iterable</i> specification provides iterator over non-unit factors only; to iterate over all factors including
 * the constant factor use {@link #iterableWithUnit()}
 *
 * @author Stanislav Poslavsky
 * @since 2.2
 */
public class Factors<E> implements Iterable<E> {
    /** The ring */
    public final Ring<E> ring;
    /** unit coefficient */
    public E unit;
    /** factors */
    public final List<E> factors;
    /** exponents */
    public final TIntArrayList exponents;

    protected Factors(Ring<E> ring, E unit, List<E> factors, TIntArrayList exponents) {
        this.ring = ring;
        this.unit = unit;
        this.factors = factors;
        this.exponents = exponents;
        if (!isUnit(unit))
            throw new IllegalArgumentException();
    }

    @Override
    public Iterator<E> iterator() {
        return factors.iterator();
    }

    /**
     * Iterator over all factors including a unit one
     *
     * @return iterator over all factors including a unit one
     */
    public Iterable<E> iterableWithUnit() {
        ArrayList<E> it = new ArrayList<>();
        if (!ring.isOne(unit))
            it.add(unit);
        it.addAll(factors);
        return it;
    }

    public boolean isUnit(E element) { return ring.isUnit(element);}

    /** Returns i-th factor */
    public E get(int i) { return factors.get(i); }

    /** Exponent of i-th factor */
    public int getExponent(int i) { return exponents.get(i); }

    /** Number of non-constant factors */
    public int size() { return factors.size(); }

    /** Whether this is a trivial factorization (contains only one factor) */
    public boolean isTrivial() { return size() == 1;}

    /** Sum all exponents */
    public int sumExponents() {
        return exponents.sum();
    }

    /** Multiply each exponent by a given factor */
    public void raiseExponents(long val) {
        for (int i = exponents.size() - 1; i >= 0; --i)
            exponents.set(i, MachineArithmetic.safeToInt(exponents.get(i) * val));
    }

    /** Sets the unit factor */
    public Factors<E> setUnit(E unit) {
        if (!isUnit(unit))
            throw new IllegalArgumentException("not a unit: " + unit);
        this.unit = unit;
        return this;
    }

    /** add another unit factor */
    public Factors<E> addUnit(E unit) {
        if (!isUnit(unit))
            throw new IllegalArgumentException("not a unit: " + unit);
        this.unit = ring.multiply(this.unit, unit);
        return this;
    }

    /** add another unit factor */
    public Factors<E> addUnit(E unit, int exponent) {
        if (!isUnit(unit))
            throw new IllegalArgumentException("not a unit: " + unit);
        this.unit = ring.multiply(this.unit, ring.pow(unit, exponent));
        return this;
    }

    /** add another factor */
    public Factors<E> addFactor(E factor, int exponent) {
        if (isUnit(factor))
            return addUnit(factor, exponent);
        factors.add(factor);
        exponents.add(exponent);
        return this;
    }

    /** add all factors from other */
    public Factors<E> addAll(Factors<E> other) {
        addUnit(other.unit);
        factors.addAll(other.factors);
        exponents.addAll(other.exponents);
        return this;
    }

    /**
     * Raise all factors to its corresponding exponents
     */
    public Factors<E> applyExponents(Ring<E> ring) {
        List<E> newFactors = new ArrayList<>();
        for (int i = 0; i < size(); i++)
            newFactors.add(ring.pow(factors.get(i), exponents.get(i)));
        return new Factors<>(ring, unit, newFactors, new TIntArrayList(ArraysUtil.arrayOf(1, size())));
    }

    /**
     * Drops constant factor from this (new instance returned)
     */
    public Factors<E> dropUnit() {
        this.unit = ring.getOne();
        return this;
    }

    /** Stream of all factors */
    public Stream<E> stream() {
        return Stream.concat(Stream.of(unit), factors.stream());
    }

    /** Stream of all factors except {@link #unit} */
    public Stream<E> streamWithoutUnit() {
        return factors.stream();
    }

    /** Array of factors without constant factor */
    public E[] toArrayWithoutUnit() {
        return factors.toArray(ring.createArray(size()));
    }

    /** Array of factors without constant factor */
    public E[] toArrayWithUnit() {
        E[] array = factors.toArray(ring.createArray(1 + size()));
        System.arraycopy(array, 0, array, 1, size());
        array[0] = unit;
        return array;
    }

    /** Multiply factors */
    public E multiply() {
        return multiply0(false);
    }

    /** Multiply with no account for exponents */
    public E multiplyIgnoreExponents() {
        return multiply0(true);
    }

    /** Square-free part */
    public E squareFreePart() {
        return multiplyIgnoreExponents();
    }

    private E multiply0(boolean ignoreExponents) {
        E r = ring.copy(unit);
        for (int i = 0; i < factors.size(); i++) {
            E tmp = ignoreExponents ? factors.get(i) : ring.pow(factors.get(i), exponents.get(i));
            r = ring.multiplyMutable(r, tmp);
        }
        return r;
    }

    /**
     * Sort factors.
     */
    public Factors<E> canonical() {
        @SuppressWarnings("unchecked")
        wrapper<E>[] wr = factors.stream().map(e -> new wrapper<>(ring, e)).toArray(wrapper[]::new);
        int[] ex = exponents.toArray();
        ArraysUtil.quickSort(wr, ex);
        factors.clear();
        exponents.clear();
        factors.addAll(Arrays.stream(wr).map(w -> w.el).collect(Collectors.toList()));
        exponents.addAll(ex);
        return this;
    }

    private static final class wrapper<E> implements Comparable<wrapper<E>> {
        final Ring<E> ring;
        final E el;

        wrapper(Ring<E> ring, E el) { this.ring = ring; this.el = el; }

        @Override
        public int compareTo(wrapper<E> o) { return ring.compare(el, o.el); }
    }

    @Override
    public String toString() {
        return toString(ToStringSupport.plain());
    }

    /**
     * @param stringifier custom stringifier
     */
    public String toString(ToStringSupport<E> stringifier) {
        return toString(stringifier, true);
    }

    /**
     * @param stringifier   custom stringifier
     * @param withExponents do stringify exponents
     */
    public String toString(ToStringSupport<E> stringifier, boolean withExponents) {
        if (factors.isEmpty())
            return stringifier.toString(unit);

        StringBuilder sb = new StringBuilder();
        if (!ring.isOne(unit)) {
            sb.append(stringifier.toString(unit));
            if (factors.size() > 0)
                sb.append("*");
        }
        for (int i = 0; ; i++) {
            sb.append("(").append(stringifier.toString(factors.get(i))).append(")");
            if (withExponents && exponents.get(i) != 1)
                sb.append("^").append(exponents.get(i));
            if (i == factors.size() - 1)
                return sb.toString();
            sb.append("*");
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Factors<?> factors1 = (Factors<?>) o;

        if (!unit.equals(factors1.unit)) return false;
        if (!factors.equals(factors1.factors)) return false;
        return exponents.equals(factors1.exponents);
    }

    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + unit.hashCode();
        result = 31 * result + factors.hashCode();
        result = 31 * result + exponents.hashCode();
        return result;
    }

    @Override
    public Factors<E> clone() {
        return new Factors<>(
                ring,
                ring.copy(unit),
                factors.stream().map(ring::copy).collect(Collectors.toList()),
                new TIntArrayList(exponents));
    }

    /** Unit factorization */
    public static <E> Factors<E> unit(Ring<E> ring, E unit) {
        if (!ring.isUnitOrZero(unit))
            throw new IllegalArgumentException("not a unit");
        return new Factors<>(ring, unit, new ArrayList<>(), new TIntArrayList());
    }

    /** Empty factorization */
    public static <E> Factors<E> empty(Ring<E> ring) {
        return unit(ring, ring.getOne());
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param ring      the ring
     * @param unit      the unit coefficient
     * @param factors   the factors
     * @param exponents the exponents
     */
    public static <E> Factors<E> of(Ring<E> ring, E unit, List<E> factors, TIntArrayList exponents) {
        if (factors.size() != exponents.size())
            throw new IllegalArgumentException();
        Factors<E> r = empty(ring).addUnit(unit);
        for (int i = 0; i < factors.size(); i++)
            r.addFactor(factors.get(i), exponents.get(i));
        return r;
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param ring    the ring
     * @param factors factors
     */
    public static <E> Factors<E> of(Ring<E> ring, E... factors) {
        return of(ring, Arrays.asList(factors));
    }

    /**
     * Factor decomposition with specified factors and exponents
     *
     * @param ring    the ring
     * @param factors factors
     */
    public static <E> Factors<E> of(Ring<E> ring, Collection<E> factors) {
        TObjectIntHashMap<E> map = new TObjectIntHashMap<>();
        for (E e : factors)
            map.adjustOrPutValue(e, 1, 1);
        List<E> l = new ArrayList<>();
        TIntArrayList e = new TIntArrayList();
        map.forEachEntry((a, b) -> {
            l.add(a);
            e.add(b);
            return true;
        });
        return of(ring, ring.getOne(), l, e);
    }
}
