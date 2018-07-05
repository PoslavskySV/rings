package cc.redberry.rings.poly;

import cc.redberry.rings.ARing;
import cc.redberry.rings.FactorDecomposition;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.io.IStringifier;
import cc.redberry.rings.poly.multivar.*;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.poly.univar.UnivariateDivision.InverseModMonomial;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Univariate quotient ring.
 *
 * @see FiniteField
 * @see AlgebraicNumberField
 */
public abstract class UnivariateQuotientRing<Poly extends IUnivariatePolynomial<Poly>>
        extends ARing<Poly> implements IPolynomialRing<Poly> {
    private static final long serialVersionUID = 1L;
    /** Irreducible polynomial that generates this finite field */
    final Poly minimalPoly;
    /** Factory polynomial */
    final Poly factory;
    /** Precomputed inverses */
    final InverseModMonomial<Poly> inverseMod;
    /** Ring cardinality */
    final BigInteger cardinality;

    /**
     * Constructs quotient ring with the specified minimal polynomial. NOTE: irreducibility test for the minimal
     * polynomial is not performed here, use {@link IrreduciblePolynomials#irreducibleQ(IUnivariatePolynomial)} to test
     * irreducibility.
     *
     * @param minimalPoly minimal polynomial
     */
    public UnivariateQuotientRing(Poly minimalPoly) {
        minimalPoly = minimalPoly.monic();
        if (minimalPoly == null)
            throw new IllegalArgumentException("Irreducible polynomial must be monic");
        this.minimalPoly = minimalPoly;
        this.factory = minimalPoly.clone();
        this.inverseMod = UnivariateDivision.fastDivisionPreConditioning(minimalPoly);
        this.cardinality = minimalPoly.coefficientRingCardinality() == null
                ? null
                : BigIntegerUtil.pow(minimalPoly.coefficientRingCardinality(), minimalPoly.degree());
    }

    /** empiric to switch between fast and plain division */
    protected boolean reduceFast(int dividendDegree) {
        int mDeg = minimalPoly.degree();
        if (dividendDegree < mDeg)
            return false;
        if (isFiniteField()) {
            if (mDeg < 8)
                return false;
            int defect = dividendDegree / mDeg;
            if (mDeg <= 20)
                return defect <= 16;
            else
                return defect <= 30;
        } else
            return false;
    }

    /**
     * Returns the irreducible polynomial that generates this finite field
     *
     * @return the irreducible polynomial that generates this finite field
     */
    public Poly getMinimalPoly() {
        return minimalPoly.clone();
    }

    public Poly element(){
        return minimalPoly.createMonomial(1);
    }
    /**
     * Gives a norm of element of this field (returned as a constant poly)
     */
    public Poly norm(Poly element) {
        return UnivariateResultants.ResultantAsPoly(minimalPoly, element);
    }

    /**
     * Gives a norm of univariate polynomial over this field
     */
    @SuppressWarnings("unchecked")
    public Poly polynomialNorm(UnivariatePolynomial<Poly> poly) {
        if (!poly.ring.equals(this))
            throw new IllegalArgumentException();
        if (minimalPoly instanceof UnivariatePolynomial)
            return (Poly) polynomialNormE((UnivariateQuotientRing) this, (UnivariatePolynomial) poly);
        if (minimalPoly instanceof UnivariatePolynomialZp64)
            return (Poly) polynomialNormZp64((UnivariateQuotientRing) this, (UnivariatePolynomial) poly);
        throw new RuntimeException();
    }

    private static <E> UnivariatePolynomial<E>
    polynomialNormE(UnivariateQuotientRing<UnivariatePolynomial<E>> ring,
                    UnivariatePolynomial<UnivariatePolynomial<E>> poly) {
        return MultivariateResultants.Resultant(
                ring.minimalPoly.asMultivariate(MonomialOrder.DEFAULT).setNVariables(2), // to bivariate
                MultivariatePolynomial.asNormalMultivariate(poly.asMultivariate(), 0),
                0).asUnivariate();
    }

    private static UnivariatePolynomialZp64
    polynomialNormZp64(UnivariateQuotientRing<UnivariatePolynomialZp64> ring,
                       UnivariatePolynomial<UnivariatePolynomialZp64> poly) {
        return MultivariateResultants.Resultant(
                ring.minimalPoly.asMultivariate(MonomialOrder.DEFAULT).setNVariables(2),
                MultivariatePolynomialZp64.asNormalMultivariate(poly.asMultivariate(), 0),
                0).asUnivariate();
    }

    /**
     * Gives a norm of multivariate polynomial over this field
     */
    @SuppressWarnings("unchecked")
    public <MPoly extends AMultivariatePolynomial>
    MPoly polynomialNorm(MultivariatePolynomial<Poly> poly) {
        if (!poly.ring.equals(this))
            throw new IllegalArgumentException();
        if (minimalPoly instanceof UnivariatePolynomial)
            return (MPoly) polynomialNormE((UnivariateQuotientRing) this, (MultivariatePolynomial) poly);
        if (minimalPoly instanceof UnivariatePolynomialZp64)
            return (MPoly) polynomialNormZp64((UnivariateQuotientRing) this, (MultivariatePolynomial) poly);
        throw new RuntimeException();
    }

    private static <E> MultivariatePolynomial<E>
    polynomialNormE(UnivariateQuotientRing<UnivariatePolynomial<E>> ring,
                    MultivariatePolynomial<UnivariatePolynomial<E>> poly) {
        return MultivariateResultants.Resultant(
                ring.minimalPoly.asMultivariate(MonomialOrder.DEFAULT).setNVariables(poly.nVariables + 1),
                MultivariatePolynomial.asNormalMultivariate(poly, 0),
                0).dropVariable(0);
    }

    private static MultivariatePolynomialZp64
    polynomialNormZp64(UnivariateQuotientRing<UnivariatePolynomialZp64> ring,
                       MultivariatePolynomial<UnivariatePolynomialZp64> poly) {
        return MultivariateResultants.Resultant(
                ring.minimalPoly.asMultivariate(MonomialOrder.DEFAULT).setNVariables(poly.nVariables + 1),
                MultivariatePolynomialZp64.asNormalMultivariate(poly, 0),
                0).dropVariable(0);
    }

    @Override
    public int nVariables() {return 1;}

    @Override
    public Poly factory() {return factory;}

    @Override
    public boolean isEuclideanRing() {return true;}

    @Override
    public BigInteger cardinality() {return cardinality;}

    @Override
    public BigInteger characteristic() {
        return minimalPoly.coefficientRingCharacteristic();
    }

    @Override
    public Poly add(Poly a, Poly b) {
        return reduceFast(Math.max(a.degree(), b.degree()))
                ? UnivariatePolynomialArithmetic.polyAddMod(a, b, minimalPoly, inverseMod, true)
                : UnivariatePolynomialArithmetic.polyAddMod(a, b, minimalPoly, true);
    }

    @Override
    public Poly subtract(Poly a, Poly b) {
        return reduceFast(Math.max(a.degree(), b.degree()))
                ? UnivariatePolynomialArithmetic.polySubtractMod(a, b, minimalPoly, inverseMod, true)
                : UnivariatePolynomialArithmetic.polySubtractMod(a, b, minimalPoly, true);
    }

    @Override
    public Poly multiply(Poly a, Poly b) {
        return reduceFast(a.degree() + b.degree())
                ? UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, minimalPoly, inverseMod, true)
                : UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, minimalPoly, true);
    }

    @Override
    public Poly negate(Poly element) {
        return reduceFast(element.degree())
                ? UnivariatePolynomialArithmetic.polyNegateMod(element, minimalPoly, inverseMod, true)
                : UnivariatePolynomialArithmetic.polyNegateMod(element, minimalPoly, true);
    }

    @Override
    public Poly addMutable(Poly a, Poly b) {
        return reduceFast(Math.max(a.degree(), b.degree()))
                ? UnivariatePolynomialArithmetic.polyAddMod(a, b, minimalPoly, inverseMod, false)
                : UnivariatePolynomialArithmetic.polyAddMod(a, b, minimalPoly, false);
    }

    @Override
    public Poly subtractMutable(Poly a, Poly b) {
        return reduceFast(Math.max(a.degree(), b.degree()))
                ? UnivariatePolynomialArithmetic.polySubtractMod(a, b, minimalPoly, inverseMod, false)
                : UnivariatePolynomialArithmetic.polySubtractMod(a, b, minimalPoly, false);
    }

    @Override
    public Poly multiplyMutable(Poly a, Poly b) {
        return reduceFast(a.degree() + b.degree())
                ? UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, minimalPoly, inverseMod, false)
                : UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, minimalPoly, false);
    }

    @Override
    public Poly negateMutable(Poly element) {
        return reduceFast(element.degree())
                ? UnivariatePolynomialArithmetic.polyNegateMod(element, minimalPoly, inverseMod, false)
                : UnivariatePolynomialArithmetic.polyNegateMod(element, minimalPoly, false);
    }

    @Override
    public Poly reciprocal(Poly element) {
        if (element.isZero())
            throw new ArithmeticException("divide by zero");
        if (isOne(element))
            return element;
        if (isMinusOne(element))
            return element;

        Poly[] xgcd = UnivariateGCD.PolynomialFirstBezoutCoefficient(element, minimalPoly);
        assert xgcd[0].isConstant();
        return xgcd[1].divideByLC(xgcd[0]);
    }

    @Override
    public FactorDecomposition<Poly> factor(Poly element) {
        return FactorDecomposition.unit(this, element);
    }

    @Override
    public Poly getZero() {
        return minimalPoly.createZero();
    }

    @Override
    public Poly getOne() {
        return minimalPoly.createOne();
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
    public Poly valueOf(long val) {
        return getOne().multiply(val);
    }

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return getOne().multiplyByBigInteger(val);
    }

    @Override
    public Poly valueOf(Poly val) {
        return reduceFast(val.degree())
                ? UnivariatePolynomialArithmetic.polyMod(val.setCoefficientRingFrom(factory), minimalPoly, inverseMod, false)
                : UnivariatePolynomialArithmetic.polyMod(val.setCoefficientRingFrom(factory), minimalPoly, false);
    }

    @Override
    public Poly copy(Poly element) {
        return element.clone();
    }

    @Override
    public Poly[] createArray(int length) {
        return minimalPoly.createArray(length);
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
        return valueOf(RandomUnivariatePolynomials.randomPoly(minimalPoly, rnd.nextInt(2 * minimalPoly.degree()), rnd));
    }

    @Override
    public Poly variable(int variable) {
        if (variable != 0)
            throw new IllegalArgumentException();
        return valueOf(minimalPoly.createMonomial(1));
    }

    @Override
    public Poly parse(String string) {
        return valueOf(factory.parsePoly(string));
    }

    /**
     * Returns iterator over all field elements
     *
     * @return iterator over all field elements
     */
    @SuppressWarnings("unchecked")
    @Override
    public Iterator<Poly> iterator() {
        if (!isFinite())
            throw new RuntimeException("Ring of infinite cardinality.");
        if (minimalPoly instanceof UnivariatePolynomial)
            return (Iterator<Poly>) new It(((UnivariatePolynomial) minimalPoly).ring, minimalPoly.degree());
        else if (minimalPoly instanceof UnivariatePolynomialZp64)
            return (Iterator<Poly>) new lIt(((UnivariatePolynomialZp64) minimalPoly).ring, minimalPoly.degree());
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        UnivariateQuotientRing<?> that = (UnivariateQuotientRing<?>) o;
        return minimalPoly.equals(that.minimalPoly);
    }

    @Override
    public int hashCode() {
        return minimalPoly.hashCode();
    }

    @Override
    public String toString(IStringifier<Poly> stringifier) {
        String cfrStr = factory.coefficientRingToString(stringifier);
        String varStr = stringifier.getBinding(factory.createMonomial(1), IStringifier.defaultVar());
        String irrStr = minimalPoly.toString(stringifier);
        return "(" + cfrStr + ")[" + varStr + "]/<" + irrStr + ">";
    }

    public String toString(String... variables) {
        return toString(IStringifier.mkPolyStringifier(factory, variables));
    }

    @Override
    public String toString() {
        return toString(IStringifier.defaultVars(1));
    }
}
