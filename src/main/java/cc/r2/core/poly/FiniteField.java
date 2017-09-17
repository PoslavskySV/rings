package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.poly.univar.*;
import cc.r2.core.poly.univar.UnivariateDivision.InverseModMonomial;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Galois field.
 *
 * @param <Poly> particular type of polynomials representing elements of this Galois field
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class FiniteField<Poly extends IUnivariatePolynomial<Poly>>
        extends ADomain<Poly> implements PolynomialDomain<Poly> {
    private static final long serialVersionUID = 1L;

    /** GF(3^3) */
    public static final FiniteField<UnivariatePolynomialZp64> GF27 = new FiniteField<>(UnivariatePolynomialZ64.create(-1, -1, 0, 1).modulus(3));
    /** GF(17^5) */
    public static final FiniteField<UnivariatePolynomialZp64> GF17p5 = new FiniteField<>(UnivariatePolynomialZ64.create(11, 11, 0, 3, 9, 9).modulus(17).monic());

    /** Irreducible polynomial that generates this finite field */
    private final Poly irreducible;
    /** Factory polynomial */
    private final Poly factory;
    /** Precomputed inverses */
    private final InverseModMonomial<Poly> inverseMod;
    /** Domain cardinality */
    private final BigInteger cardinality;

    /**
     * Constructs finite field from the specified irreducible polynomial. NOTE: irreducibility test for the input
     * polynomial is not performed here, use {@link cc.r2.core.poly.univar.IrreduciblePolynomials#irreducibleQ(IUnivariatePolynomial)}
     * to test irreducibility.
     *
     * @param irreducible irreducible polynomial
     */
    public FiniteField(Poly irreducible) {
        if (!irreducible.isOverField())
            throw new IllegalArgumentException();
        this.irreducible = irreducible;
        this.factory = irreducible.clone();
        this.inverseMod = UnivariateDivision.fastDivisionPreConditioning(irreducible);
        this.cardinality = BigIntegerArithmetics.pow(irreducible.coefficientDomainCardinality(), irreducible.degree());
    }

    /**
     * Returns the irreducible polynomial that generates this finite field
     *
     * @return the irreducible polynomial that generates this finite field
     */
    public Poly getIrreducible() {
        return irreducible.clone();
    }

    @Override
    public int nVariables() { return 1; }

    @Override
    public Poly factory() { return factory; }

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public BigInteger cardinality() {
        return cardinality;
    }

    @Override
    public BigInteger characteristic() {
        return irreducible.coefficientDomainCharacteristics();
    }

    @Override
    public Poly add(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polyAddMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly subtract(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polySubtractMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly multiply(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly negate(Poly element) {
        return UnivariatePolynomialArithmetic.polyNegateMod(element, irreducible, inverseMod, true);
    }

    @Override
    public Poly addMutable(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polyAddMod(a, b, irreducible, inverseMod, false);
    }

    @Override
    public Poly subtractMutable(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polySubtractMod(a, b, irreducible, inverseMod, false);
    }

    @Override
    public Poly multiplyMutable(Poly a, Poly b) {
        return UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, irreducible, inverseMod, false);
    }

    @Override
    public Poly negateMutable(Poly element) {
        return UnivariatePolynomialArithmetic.polyNegateMod(element, irreducible, inverseMod, false);
    }

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {
        return a.createArray(multiply(a, reciprocal(b)), getZero());
    }

    @Override
    public Poly remainder(Poly dividend, Poly divider) {
        return getZero();
    }

    @Override
    public Poly reciprocal(Poly element) {
        if (element.isZero())
            throw new ArithmeticException("divide by zero");
        Poly
                t = getZero(),
                newt = getOne(),
                r = irreducible.clone(),
                newr = element.clone();
        Poly tmp;
        while (!newr.isZero()) {
            Poly quotient = UnivariateDivision.quotient(r, newr, true);

            tmp = r;
            r = newr;
            newr = tmp.clone().subtract(quotient.clone().multiply(newr));

            tmp = t;
            t = newt;
            newt = tmp.clone().subtract(quotient.clone().multiply(newt));
        }
        if (r.degree() > 0)
            throw new IllegalArgumentException("Either p is not irreducible or a is a multiple of p, p =  " + irreducible + ", a = " + element);
        return t.divideByLC(r);
    }

    @Override
    public Poly gcd(Poly a, Poly b) {
        return UnivariateGCD.PolynomialGCD(a, b);
    }

    @Override
    public Poly getZero() {
        return irreducible.createZero();
    }

    @Override
    public Poly getOne() {
        return irreducible.createOne();
    }

    @Override
    public boolean isZero(Poly element) {
        return element.isZero();
    }

    @Override
    public boolean isOne(Poly element) {
        return element.isOne();
    }

    @Override
    public boolean isUnit(Poly element) {
        return !isZero(element);
    }

    @Override
    public Poly valueOf(long val) {
        return getOne().multiply(val);
    }

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return getOne().multiplyByBigInteger(val);
    }

    @Override
    public Poly valueOf(Poly val) {
        return UnivariatePolynomialArithmetic.polyMod(val, irreducible, inverseMod, true);
    }

    @Override
    public Poly copy(Poly element) {
        return element.clone();
    }

    @Override
    public Poly[] createArray(int length) {
        return irreducible.createArray(length);
    }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[][] createArray2d(int length) {
        Poly[] array = createArray(0);
        return (Poly[][]) Array.newInstance(array.getClass(), length);
    }

    @Override
    public Poly[][] createArray2d(int m, int n) {
        Poly[][] arr = createArray2d(m);
        for (int i = 0; i < arr.length; i++)
            arr[i] = createArray(n);
        return arr;
    }

    @Override
    public int compare(Poly o1, Poly o2) {
        return o1.compareTo(o2);
    }

    @Override
    public Poly randomElement(RandomGenerator rnd) {
        return valueOf(RandomUnivariatePolynomials.randomPoly(irreducible, rnd.nextInt(2 * irreducible.degree()), rnd));
    }

    @Override
    public Poly variable(int variable) {
        if (variable != 0)
            throw new IllegalArgumentException();
        return valueOf(irreducible.createMonomial(1));
    }

    @Override
    public Poly parse(String string) {
        return valueOf(irreducible.parsePoly(string));
    }

    @Override
    public Poly parse(String string, String[] variables) {
        return valueOf(irreducible.parsePoly(string, variables));
    }

    /**
     * Returns iterator over all field elements
     *
     * @return iterator over all field elements
     */
    @SuppressWarnings("unchecked")
    @Override
    public Iterator<Poly> iterator() {
        if (irreducible instanceof UnivariatePolynomial)
            return (Iterator<Poly>) new It(((UnivariatePolynomial) irreducible).domain, irreducible.degree());
        else if (irreducible instanceof UnivariatePolynomialZp64)
            return (Iterator<Poly>) new lIt(((UnivariatePolynomialZp64) irreducible).domain, irreducible.degree());
        throw new RuntimeException();
    }

    private static final class It<E> implements Iterator<UnivariatePolynomial<E>> {
        final Domain<E> domain;
        final E[] data;
        final Iterator<E>[] iterators;

        @SuppressWarnings("unchecked")
        It(Domain<E> domain, int degree) {
            this.domain = domain;
            this.data = domain.createArray(degree);
            this.iterators = new Iterator[degree];
            for (int i = 0; i < iterators.length; i++)
                iterators[i] = domain.iterator();
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
                return UnivariatePolynomial.create(domain, data.clone());
            }
            int i = 0;
            if (!iterators[i].hasNext())
                while (i < iterators.length && !iterators[i].hasNext()) {
                    iterators[i] = domain.iterator();
                    data[i] = iterators[i].next();
                    ++i;
                }

            if (i >= iterators.length)
                return null;

            data[i] = iterators[i].next();
            return UnivariatePolynomial.create(domain, data.clone());
        }
    }

    private static final class lIt implements Iterator<UnivariatePolynomialZp64> {
        final IntegersZp64 domain;
        final long[] data;

        @SuppressWarnings("unchecked")
        lIt(IntegersZp64 domain, int degree) {
            this.domain = domain;
            this.data = new long[degree];
        }

        @Override
        public boolean hasNext() {
            return Arrays.stream(data).anyMatch(l -> l < (domain.modulus - 1));
        }

        private boolean first = true;

        @Override
        public UnivariatePolynomialZp64 next() {
            if (first) {
                first = false;
                return UnivariatePolynomialZp64.createUnsafe(domain, data.clone());
            }
            int i = 0;
            if (data[i] >= domain.modulus - 1)
                while (i < data.length && data[i] >= domain.modulus - 1) {
                    data[i] = 0;
                    ++i;
                }

            if (i >= data.length)
                return null;

            ++data[i];
            return UnivariatePolynomialZp64.createUnsafe(domain, data.clone());
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FiniteField<?> that = (FiniteField<?>) o;

        if (!irreducible.equals(that.irreducible)) return false; // todo: need to check this?
        return cardinality.equals(that.cardinality);
    }

    @Override
    public int hashCode() {
        int result = irreducible.hashCode();
        result = 31 * result + cardinality.hashCode();
        return result;
    }

    @Override
    public String toString(String[] variables) {
        return toString(irreducible.coefficientDomainToString(), variables);
    }

    @Override
    public String toString(String coefficientDomain, String[] variables) {
        return "(" + coefficientDomain + ")[" + variables[0] + "]/<" + irreducible.toString(variables) + ">";
    }

    public String toString(String coefficientDomain, ToStringSupport<Poly> irreducibleToString, String[] variables) {
        return "(" + coefficientDomain + ")[" + variables[0] + "]/<" + irreducibleToString.toString(irreducible) + ">";
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(1));
    }
}
