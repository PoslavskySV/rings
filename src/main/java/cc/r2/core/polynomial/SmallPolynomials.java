package cc.r2.core.polynomial;


import cc.r2.core.number.primes.PrimesIterator;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static cc.r2.core.polynomial.LongArithmetics.*;
import static cc.r2.core.polynomial.SmallPolynomialArithmetics.derivative;

public final class SmallPolynomials {
    private SmallPolynomials() {
    }

    /**
     * Plain Euclidean algorithm which fails if intermediate polynomials are not divisible at some step
     *
     * @param a poly
     * @param b poly
     * @return a list of polynomial remainders where the last element is GCD
     * @throws IllegalArgumentException if at some step intermediate polynomials are not divisible
     */
    public static PolynomialRemainders Euclid(final MutableLongPoly a,
                                              final MutableLongPoly b) {
        if (a.degree < b.degree)
            return Euclid(b, a);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(a.clone()); prs.add(b.clone());

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(prs);

        MutableLongPoly x = a, y = b, r;
        while (true) {
            //TODO replace with remainder!!!
            MutableLongPoly[] tmp = divideAndRemainder(x, y, true);
            if (tmp == null)
                throw new IllegalArgumentException("Not divisible: (" + x + ") / (" + y + ")");

            r = tmp[1];
            if (r.isZero())
                break;

            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders(prs);
    }

    /**
     * Euclidean algorithm
     *
     * @param a poly
     * @param b poly
     * @return a list of polynomial remainders where the last element is GCD
     */
    public static PolynomialRemainders Euclid(final MutableLongPoly a,
                                              final MutableLongPoly b,
                                              final long modulus) {
        if (a.degree < b.degree)
            return Euclid(b, a, modulus);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(a.clone().modulus(modulus)); prs.add(b.clone().modulus(modulus));

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(prs);

        MutableLongPoly x = a, y = b, r;
        while (true) {
            //TODO replace with remainder!!!
            MutableLongPoly[] tmp = divideAndRemainder(x, y, modulus, true);
            assert tmp != null;
            r = tmp[1];
            if (r.isZero())
                break;
            prs.add(r);
            x = y;
            y = r;
        }
        return new PolynomialRemainders(prs);
    }

    /**
     * Euclidean algorithm for polynomials that uses pseudo division
     *
     * @param a            poly
     * @param b            poly
     * @param primitivePRS whether to use primitive polynomial remainders
     * @return a list of polynomial remainders where the last element is GCD
     */
    public static PolynomialRemainders PolynomialEuclid(final MutableLongPoly a,
                                                        final MutableLongPoly b,
                                                        boolean primitivePRS) {
        if (a.degree < b.degree)
            return PolynomialEuclid(b, a, primitivePRS);


        if (a.isZero() || b.isZero()) return new PolynomialRemainders(a.clone(), b.clone());

        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        MutableLongPoly aPP = a.clone().divide(aContent), bPP = b.clone().divide(bContent);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        MutableLongPoly x = aPP, y = bPP, r;
        while (true) {
            //TODO replace with remainder!!!
            MutableLongPoly[] tmp = pseudoDivideAndRemainder(x, y, true);
            assert tmp != null;
            r = tmp[1];
            if (r.isZero())
                break;
            if (primitivePRS)
                r = primitivePart(r);
            prs.add(r);
            x = y;
            y = r;
        }
        PolynomialRemainders res = new PolynomialRemainders(prs);
        primitivePart(res.gcd()).multiply(contentGCD);
        return res;
    }

    /**
     * Euclidean algorithm for polynomials which produces subresultants sequence
     *
     * @param a poly
     * @param b poly
     * @return subresultant sequence where the last element is GCD
     */
    public static PolynomialRemainders SubresultantEuclid(final MutableLongPoly a,
                                                          final MutableLongPoly b) {
        if (b.degree > a.degree)
            return SubresultantEuclid(b, a);

        if (a.isZero() || b.isZero()) return new PolynomialRemainders(a.clone(), b.clone());


        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        MutableLongPoly aPP = a.clone().divide(aContent), bPP = b.clone().divide(bContent);

        ArrayList<MutableLongPoly> prs = new ArrayList<>();
        prs.add(aPP); prs.add(bPP);

        TLongArrayList beta = new TLongArrayList(), psi = new TLongArrayList();
        TIntArrayList deltas = new TIntArrayList();

        long cBeta, cPsi;
        for (int i = 0; ; i++) {
            MutableLongPoly curr = prs.get(i);
            MutableLongPoly next = prs.get(i + 1);
            int delta = curr.degree - next.degree;
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? 1 : -1;
                cPsi = -1;
            } else {
                cPsi = pow(-curr.lc(), deltas.get(i - 1));
                if (deltas.get(i - 1) < 1) {
                    cPsi = multiply(cPsi, pow(psi.get(i - 1), -deltas.get(i - 1) + 1));
                } else {
                    long tmp = pow(psi.get(i - 1), deltas.get(i - 1) - 1);
                    assert cPsi % tmp == 0;
                    cPsi /= tmp;
                }
                cBeta = multiply(-curr.lc(), pow(cPsi, delta));
            }

            //TODO replace with remainder!!!
            MutableLongPoly q = pseudoDivideAndRemainder(curr, next, true)[1];
            if (q.isZero())
                break;

            q = q.divide(cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        PolynomialRemainders res = new PolynomialRemainders(prs);
        res.gcd().multiply(contentGCD);
        return res;
    }


    /**
     * Representation for polynomial remainders sequence produced by the Euclidean algorithm
     */
    public static final class PolynomialRemainders {
        public final ArrayList<MutableLongPoly> remainders;

        public PolynomialRemainders(MutableLongPoly... remainders) {
            this(new ArrayList<>(Arrays.asList(remainders)));
        }

        public PolynomialRemainders(ArrayList<MutableLongPoly> remainders) {
            this.remainders = remainders;
        }

        public MutableLongPoly gcd() {
            if (remainders.size() == 2 && remainders.get(1).isZero())
                return remainders.get(0);
            return remainders.get(remainders.size() - 1);
        }
    }

    /**
     * Modular GCD algorithm for polynomials
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static MutableLongPoly ModularGCD(MutableLongPoly a, MutableLongPoly b) {
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree < b.degree)
            return ModularGCD(b, a);
        long aContent = content(a), bContent = content(b);
        long contentGCD = gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return MutableLongPoly.create(contentGCD);

        return ModularGCD0(a.clone().divide(aContent), b.clone().divide(bContent)).multiply(contentGCD);

    }

    /** modular GCD for primitive polynomials */
    @SuppressWarnings("ConstantConditions")
    private static MutableLongPoly ModularGCD0(MutableLongPoly a, MutableLongPoly b) {
        if (a.degree < b.degree)
            return ModularGCD(b, a);

        long lcGCD = gcd(a.lc(), b.lc());
        double bound = Math.sqrt(a.degree + 1) * (1L << a.degree) * Math.max(a.norm(), b.norm()) * lcGCD;

        MutableLongPoly previousBase, base = null;
        long basePrime = -1;

        PrimesIterator primesLoop = new PrimesIterator(3);
        while (true) {
            long prime = primesLoop.take();
            assert prime != -1 : "long overflow";

            if (a.lc() % prime == 0 || b.lc() % prime == 0)
                continue;

            MutableLongPoly modularGCD = Euclid(a, b, prime).gcd();
            //clone if necessary
            if (modularGCD == a || modularGCD == b)
                modularGCD = modularGCD.clone();

            //coprime polynomials
            if (modularGCD.degree == 0)
                return MutableLongPoly.one();


            //save the base
            if (base == null) {
                //make base monic and multiply lcGCD
                modularGCD.monic(lcGCD, prime);
                base = modularGCD;
                basePrime = prime;
                continue;
            }

            //unlucky base => start over
            if (base.degree > modularGCD.degree) {
                base = null;
                basePrime = -1;
                continue;
            }

            //skip unlucky prime
            if (base.degree < modularGCD.degree)
                continue;

            //cache current base
            previousBase = base.clone();

            //lifting
            long newBasePrime = multiply(basePrime, prime);
            long monicFactor = modInverse(modularGCD.lc(), prime);
            long lcMod = mod(lcGCD, prime);
            for (int i = 0; i <= base.degree; ++i) {
                //this is monic modularGCD multiplied by lcGCD mod prime
                long oth = mod(multiply(mod(multiply(modularGCD.data[i], monicFactor), prime), lcMod), prime);
                base.data[i] = ChineseRemainders(basePrime, prime, base.data[i], oth);
            }
            base.fixDegree();
            basePrime = newBasePrime;

            //either trigger Mignotte's bound or two trials didn't change the result, probably we are done
            if ((double) basePrime >= 2 * bound || base.equals(previousBase)) {
                MutableLongPoly candidate = primitivePart(base.clone().symModulus(basePrime));
                //first check b since b is less degree
                MutableLongPoly[] div;
                div = divideAndRemainder(b, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                div = divideAndRemainder(a, candidate, true);
                if (div == null || !div[1].isZero())
                    continue;

                return candidate;
            }
        }
    }

    /**
     * Computes GCD of two polynomials
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    public static MutableLongPoly PolynomialGCD(MutableLongPoly a, MutableLongPoly b) {
        return ModularGCD(a, b);
    }

    /**
     * Computes GCD of two polynomials modulo prime
     *
     * @param a       the first polynomial
     * @param b       the second polynomial
     * @param modulus prime modulus
     * @return GCD of two polynomials
     */
    public static MutableLongPoly PolynomialGCD(MutableLongPoly a, MutableLongPoly b, long modulus) {
        return Euclid(a, b, modulus).gcd();
    }

    /**
     * Performs square-free factorization of a poly.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static Factorization SquareFreeFactorization(MutableLongPoly poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun(poly);

        MutableLongPoly expFree = MutableLongPoly.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
        return SquareFreeFactorizationYun(expFree).addFactor(MutableLongPoly.createMonomial(1, exponent), 1);
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static Factorization SquareFreeFactorizationYun(MutableLongPoly poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = content(poly);
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divide(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutableLongPoly derivative = derivative(poly), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutableLongPoly
                quot = divideAndRemainder(poly, gcd, false)[0], // safely destroy (cloned) poly (not used further)
                dQuot = divideAndRemainder(derivative, gcd, false)[0]; // safely destroy (cloned) derivative (not used further)

        ArrayList<MutableLongPoly> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            dQuot = dQuot.subtract(derivative(quot));
            MutableLongPoly factor = PolynomialGCD(quot, dQuot);
            quot = divideAndRemainder(quot, factor, false)[0]; // can destroy quot in divideAndRemainder
            dQuot = divideAndRemainder(dQuot, factor, false)[0]; // can destroy dQuot in divideAndRemainder
            if (!factor.isOne()) {
                factors.add(factor);
                exponents.add(i);
            }
        }
        return new Factorization(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a poly using Musser's algorithm
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static Factorization SquareFreeFactorizationMusser(MutableLongPoly poly) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long content = content(poly);
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divide(content);
        if (poly.degree <= 1)
            return oneFactor(poly, content);

        MutableLongPoly derivative = derivative(poly), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant())
            return oneFactor(poly, content);

        MutableLongPoly quot = divideAndRemainder(poly, gcd, false)[0]; // safely destroy (cloned) poly

        ArrayList<MutableLongPoly> factors = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        int i = 0;
        while (true) {
            ++i;
            MutableLongPoly nextQuot = PolynomialGCD(gcd, quot);
            gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // safely destroy gcd (reassigned)
            MutableLongPoly factor = divideAndRemainder(quot, nextQuot, false)[0]; // safely destroy quot (reassigned further)
            if (!factor.isConstant()) {
                factors.add(factor);
                exponents.add(i);
            }
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
        return new Factorization(factors, exponents, content);
    }

    /**
     * Performs square-free factorization of a {@code poly} modulo {@code modulus}.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return square-free decomposition modulo {@code modulus}
     */
    public static Factorization SquareFreeFactorization(MutableLongPoly poly, long modulus) {
        if (modulus >= Integer.MAX_VALUE)// <- just to be on the safe side)
            throw new IllegalArgumentException();

        poly = poly.clone().modulus(modulus);
        long lc = poly.lc();
        //make poly monic
        poly = poly.monic(modulus);

        if (poly.isConstant())
            return oneFactor(lc);

        if (poly.degree <= 1)
            return oneFactor(poly, lc);

        Factorization factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly, modulus);
        else {
            MutableLongPoly expFree = MutableLongPoly.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
            factorization = SquareFreeFactorizationMusser0(expFree, modulus).addFactor(MutableLongPoly.createMonomial(1, exponent), 1);
        }

        return factorization.setFactor(lc);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static Factorization SquareFreeFactorizationMusser0(MutableLongPoly poly, long modulus) {
        poly.monic(modulus);
        if (poly.isConstant())
            return oneFactor(poly.lc());

        if (poly.degree <= 1)
            return oneFactor(poly, 1);

        MutableLongPoly derivative = derivative(poly, modulus);
        if (!derivative.isZero()) {
            MutableLongPoly gcd = PolynomialGCD(poly, derivative, modulus);
            if (gcd.isConstant())
                return oneFactor(poly, 1);
            MutableLongPoly quot = divideAndRemainder(poly, gcd, modulus, false)[0]; // can safely destroy poly (not used further)

            ArrayList<MutableLongPoly> factors = new ArrayList<>();
            TIntArrayList exponents = new TIntArrayList();
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                MutableLongPoly nextQuot = PolynomialGCD(gcd, quot, modulus);
                MutableLongPoly factor = divideAndRemainder(quot, nextQuot, modulus, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant()) {
                    factors.add(factor.monic(modulus));
                    exponents.add(i);
                }
                gcd = divideAndRemainder(gcd, nextQuot, modulus, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd, modulus);
                Factorization gcdFactorization = SquareFreeFactorizationMusser0(gcd, modulus);
                gcdFactorization.raiseExponents((int) modulus);
                MutableLongPoly[] newFactors = factors.toArray(new MutableLongPoly[gcdFactorization.factors.length + factors.size()]);
                int[] newExponents = Arrays.copyOf(exponents.toArray(), gcdFactorization.factors.length + factors.size());
                System.arraycopy(gcdFactorization.factors, 0, newFactors, factors.size(), gcdFactorization.factors.length);
                System.arraycopy(gcdFactorization.exponents, 0, newExponents, exponents.size(), gcdFactorization.exponents.length);
                return new Factorization(newFactors, newExponents, 1);
            } else
                return new Factorization(factors.toArray(new MutableLongPoly[factors.size()]), exponents.toArray(), 1);
        } else {
            MutableLongPoly pRoot = pRoot(poly, modulus);
            Factorization factorization = SquareFreeFactorizationMusser0(pRoot, modulus);
            factorization.raiseExponents((int) modulus);
            return new Factorization(factorization.factors, factorization.exponents, 1);
        }
    }

    /** p-th root of poly */
    private static MutableLongPoly pRoot(MutableLongPoly poly, long modulus) {
        assert poly.degree % modulus == 0;
        long[] rootData = new long[poly.degree / (int) modulus + 1];
        Arrays.fill(rootData, 0);
        for (int i = poly.degree; i >= 0; --i)
            if (poly.data[i] != 0) {
                assert i % modulus == 0;
                rootData[i / (int) modulus] = poly.data[i];
            }
        return MutableLongPoly.create(rootData);
    }

    /**
     * Returns {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if and only if {@code poly} is square-free and {@code false} otherwise
     */
    public static boolean isSquareFree(MutableLongPoly poly) {
        return PolynomialGCD(poly, derivative(poly)).isConstant();
    }

    /**
     * Returns {@code true} if and only if {@code poly} is square-free modulo {@code modulus} and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if and only if {@code poly} is square-free modulo {@code modulus} and {@code false} otherwise
     */
    public static boolean isSquareFree(MutableLongPoly poly, long modulus) {
        return PolynomialGCD(poly, derivative(poly, modulus), modulus).isConstant();
    }

//    /** Berlekamp's Q-matrix */
//    public static final class QMatrix {
//        final long[] data;
//        final int size;
//
//        public QMatrix(int degree) {
//            this.size = degree - 1;
//            this.data = new long[size * size];
//        }
//
//        void set(int col, int row, long val) {
//            data[row + size * col] = val;
//        }
//
//        long get(int col, int row) {
//            return data[row + size * col];
//        }
//    }
//
//
//    static long[][] QMatrix(MutableLongPoly poly, long modulus) {
//        int pDegree = poly.degree;
//        long[][] matrix = new long[pDegree][pDegree];
//        long[] prevRow = new long[pDegree];
//        prevRow[0] = 1;
//        matrix[0] = prevRow;
//        for (int i = 1; i <= (pDegree - 1) * modulus; i++) {
//            long nextRow[] = new long[pDegree];
//            nextRow[0] = symMod(-prevRow[pDegree - 1] * poly.data[0], modulus);
//            for (int j = 1; j < poly.degree; j++) {
//                nextRow[j] = symMod(prevRow[j - 1] - prevRow[pDegree - 1] * poly.data[j], modulus);
//            }
//            if (i % modulus == 0)
//                matrix[i / (int) modulus] = nextRow.clone();
//            prevRow = nextRow;
//        }
//        return matrix;
//    }

    /**
     * Performs distinct-degree factorization for square-free polynomial {@code poly} modulo {@code modulus}.
     * In the case of not square-free input, the algorithm works, but the resulting factorization is not compete.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    public static Factorization DistinctDegreeFactorization(MutableLongPoly poly, long modulus) {
        if (poly.isConstant())
            return oneFactor(poly.lc());

        long factor = mod(poly.lc(), modulus);
        MutableLongPoly base = poly.clone().monic(modulus);
        MutableLongPoly polyModulus = base.clone();

        if (base.degree <= 1)
            return oneFactor(base, factor);

        if (base.isMonomial())
            return oneFactor(base, factor);

        MutableLongPoly exponent = MutableLongPoly.create(0, 1);
        ArrayList<MutableLongPoly> factors = new ArrayList<>();
        TIntArrayList degrees = new TIntArrayList();
        int i = 0;
        long total = 0;
        while (!base.isConstant()) {
//            System.out.println(i);
            ++i;
//            System.out.println("before powmod");
//            System.out.println("Exponent:" + exponent.toStringForCopy());
//            System.out.println("Modulus:" + modulus);
//            System.out.println("PolyModulus:" + polyModulus.toStringForCopy());
            long start = System.nanoTime();
            exponent = SmallPolynomialArithmetics.powMod(exponent, modulus, base, modulus, false);
            total += (System.nanoTime() - start);
//            System.out.println("after powmod");
            MutableLongPoly tmpExponent = exponent.clone();
            //subtract x^1 effectively
            tmpExponent.ensureCapacity(1);
            tmpExponent.data[1] = subtractMod(tmpExponent.data[1], 1, modulus);
            tmpExponent.fixDegree();
//            System.out.println("before gcd");
            MutableLongPoly gcd = PolynomialGCD(tmpExponent, base, modulus);
//            System.out.println("after gcd");
            if (!gcd.isConstant()) {
                factors.add(gcd.monic(modulus));
                degrees.add(i);
            }
            assert divideAndRemainder(base, gcd, modulus, true)[1].isZero();
            base = divideAndRemainder(base, gcd, modulus, false)[0]; //can safely destroy reused base

            if (base.degree < 2 * (i + 1)) {// <- early termination
                if (!base.isConstant()) {
                    factors.add(base.monic(modulus));
                    degrees.add(i + 1);
                }
                break;
            }
        }

//        System.out.println(total);
        return new Factorization(factors, degrees, factor);
    }

    /**
     * Performs square-free factorization followed by distinct-degree factorization modulo {@code modulus}.
     *
     * @param poly    the polynomial
     * @param modulus prime modulus
     * @return square-free and distinct-degree decomposition of {@code poly} modulo {@code modulus}
     */
    @SuppressWarnings("ConstantConditions")
    public static Factorization DistinctDegreeFactorizationComplete(MutableLongPoly poly, long modulus) {
        Factorization squareFree = SquareFreeFactorization(poly, modulus);
        long overallFactor = squareFree.factor;
        ArrayList<MutableLongPoly> finalFactors = new ArrayList<>(squareFree.factors.length);
        TIntArrayList finalExponents = new TIntArrayList(squareFree.factors.length);
        for (int i = squareFree.factors.length - 1; i >= 0; --i) {
            Factorization dd = DistinctDegreeFactorization(squareFree.factors[i], modulus);
            int nFactors = dd.factors.length;
            finalFactors.ensureCapacity(finalFactors.size() + nFactors);
            finalExponents.ensureCapacity(finalExponents.size() + nFactors);
            for (int j = nFactors - 1; j >= 0; --j) {
                finalFactors.add(dd.factors[j]);
                finalExponents.add(squareFree.exponents[i]);
            }
            overallFactor = multiplyMod(overallFactor, dd.factor, modulus);
        }

        return new Factorization(finalFactors, finalExponents, overallFactor);
    }

    /** pretty print for factorization */
    private static String toStringFactorization(Object[] factors, int[] exponents, long factor, boolean infoAsExponents) {
        StringBuilder sb = new StringBuilder();
        if (factor != 1) {
            sb.append(factor);
            if (factors.length > 0)
                sb.append("*");
        }
        for (int i = 0; ; i++) {
            sb.append("(").append(factors[i]).append(")");
            if (infoAsExponents && exponents[i] != 1)
                sb.append("^").append(exponents[i]);
            if (i == factors.length - 1)
                return sb.toString();
            sb.append("*");
        }
    }

    static Factorization oneFactor(long factor) {
        return new Factorization(new MutableLongPoly[0], new int[0], factor);
    }

    static Factorization oneFactor(MutableLongPoly poly, long factor) {
        return new Factorization(new MutableLongPoly[]{poly}, new int[]{1}, factor);
    }

    /**
     * Polynomial factorization
     */
    public static final class Factorization {
        /** Integer factor (polynomial content) */
        final long factor;
        /** Factors */
        final MutableLongPoly[] factors;
        /** Either exponents or distinct-degree powers */
        final int[] exponents;

        Factorization(MutableLongPoly[] factors, int[] exponents, long factor) {
            this.factors = factors;
            this.exponents = exponents;
            this.factor = factor;
        }

        Factorization(List<MutableLongPoly> factors, TIntArrayList exponents, long factor) {
            this(factors.toArray(new MutableLongPoly[factors.size()]), exponents.toArray(), factor);
        }

        protected String toString(boolean infoAsExponents) {
            return toStringFactorization(factors, exponents, factor, infoAsExponents);
        }

        Factorization setFactor(long factor) {
            if (factor == this.factor) return this;
            return new Factorization(factors, exponents, factor);
        }

        void raiseExponents(int val) {
            for (int i = exponents.length - 1; i >= 0; --i)
                exponents[i] *= val;
        }

        Factorization canonical() {
            MutableLongPoly[] fTmp = factors.clone();
            int[] eTmp = exponents.clone();
            for (int i = fTmp.length - 1; i >= 0; --i) {
                MutableLongPoly poly = fTmp[i];
                if (poly.isMonomial() && eTmp[i] != 1) {
                    int degree = poly.degree;
                    poly.ensureCapacity(poly.degree * eTmp[i]);
                    poly.data[degree * eTmp[i]] = poly.data[degree];
                    poly.data[degree] = 0;
                    eTmp[i] = 1;
                    assert poly.isMonomial();
                }
            }

            ArraysUtil.quickSort(fTmp, eTmp);
            return new Factorization(fTmp, eTmp, factor);
        }

        Factorization addFactor(MutableLongPoly poly, int exponent) {
            MutableLongPoly[] factors = Arrays.copyOf(this.factors, this.factors.length + 1);
            int[] exponents = Arrays.copyOf(this.exponents, this.exponents.length + 1);
            factors[factors.length - 1] = poly;
            exponents[exponents.length - 1] = exponent;
            return new Factorization(factors, exponents, factor);
        }

        @Override
        public String toString() {
            return toString(true);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Factorization that = (Factorization) o;

            if (factor != that.factor) return false;
            // Probably incorrect - comparing Object[] arrays with Arrays.equals
            if (!Arrays.equals(factors, that.factors)) return false;
            return Arrays.equals(exponents, that.exponents);
        }

        @Override
        public int hashCode() {
            int result = (int) (factor^(factor >>> 32));
            result = 31 * result + Arrays.hashCode(factors);
            result = 31 * result + Arrays.hashCode(exponents);
            return result;
        }
    }

    /**
     * Returns the content of the poly
     *
     * @param poly poly
     * @return polynomial content
     */
    public static long content(MutableLongPoly poly) {
        if (poly.degree == 0)
            return poly.data[0];
        return gcd(poly.data, 0, poly.degree + 1);
    }

    /**
     * Reduces poly to its primitive part
     *
     * @param poly polynomial
     * @return primitive part (poly will be modified)
     */
    public static MutableLongPoly primitivePart(MutableLongPoly poly) {
        long content = content(poly);
        if (content == 1)
            return poly;
        if (poly.lc() < 0)
            content = -content;
        for (int i = 0; i <= poly.degree; ++i) {
            assert poly.data[i] % content == 0;
            poly.data[i] = poly.data[i] / content;
        }
        return poly;
    }

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       boolean copy) {
        if (dividend.isZero())
            return new MutableLongPoly[]{MutableLongPoly.zero(), MutableLongPoly.zero()};
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0) {
            MutableLongPoly div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new MutableLongPoly[]{div, MutableLongPoly.zero()};
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
        return divideAndRemainderGeneral0(dividend, divider, 1, copy);
    }

    /**
     * Returns quotient and remainder using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] pseudoDivideAndRemainder(final MutableLongPoly dividend,
                                                             final MutableLongPoly divider,
                                                             final boolean copy) {
        if (dividend.isZero())
            return new MutableLongPoly[]{MutableLongPoly.zero(), MutableLongPoly.zero()};
        if (dividend.degree < divider.degree)
            return null;
        long factor = pow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new MutableLongPoly[]{(copy ? dividend.clone() : dividend).multiply(factor / dividend.lc()), MutableLongPoly.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderGeneral0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static MutableLongPoly[] divideAndRemainderGeneral0(final MutableLongPoly dividend,
                                                        final MutableLongPoly divider,
                                                        final long dividendRaiseFactor,
                                                        final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor),
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0)
                    return null;

                quotient.data[i] = remainder.lc() / divider.lc();
                remainder.subtract(divider, quotient.data[i], i);

            } else quotient.data[i] = 0;
        }

        quotient.fixDegree();
        return new MutableLongPoly[]{quotient, remainder};
    }

    /**
     * Returns quotient and remainder using adaptive pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    static MutableLongPoly[] pseudoDivideAndRemainderAdaptive(final MutableLongPoly dividend,
                                                              final MutableLongPoly divider,
                                                              final boolean copy) {
        if (dividend.isZero())
            return new MutableLongPoly[]{MutableLongPoly.zero(), MutableLongPoly.zero()};
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0)
            return new MutableLongPoly[]{copy ? dividend.clone() : dividend, MutableLongPoly.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static MutableLongPoly[] pseudoDivideAndRemainderAdaptive0(final MutableLongPoly dividend,
                                                               final MutableLongPoly divider,
                                                               final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly
                remainder = copy ? dividend.clone() : dividend,
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0) {
                    long gcd = gcd(remainder.lc(), divider.lc());
                    long factor = divider.lc() / gcd;
                    remainder.multiply(factor);
                    quotient.multiply(factor);
                }

                quotient.data[i] = remainder.lc() / divider.lc();
                remainder.subtract(divider, quotient.data[i], i);

            } else quotient.data[i] = 0;
        }

        quotient.fixDegree();
        return new MutableLongPoly[]{quotient, remainder};
    }

    /**
     * Returns quotient and remainder modulo {@code modulus}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutableLongPoly[] divideAndRemainder(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       final long modulus,
                                                       final boolean copy) {
        if (dividend.isZero())
            return new MutableLongPoly[]{MutableLongPoly.zero(), MutableLongPoly.zero()};
        if (dividend.degree < divider.degree)
            return null;
        if (divider.degree == 0)
            return new MutableLongPoly[]{(copy ? dividend.clone() : dividend).multiply(modInverse(divider.lc(), modulus), modulus), MutableLongPoly.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, modulus, copy);
        return divideAndRemainderModulus(dividend, divider, modulus, copy);
    }

    /** Plain school implementation */
    static MutableLongPoly[] divideAndRemainderModulus(final MutableLongPoly dividend,
                                                       final MutableLongPoly divider,
                                                       final long modulus,
                                                       final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly
                remainder = copy ? dividend.clone() : dividend,
                quotient = new MutableLongPoly(dividend.degree - divider.degree);

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient.data[i] = divideMod(remainder.lc(), divider.lc(), modulus);
                remainder.subtract(divider, quotient.data[i], i, modulus);
            } else quotient.data[i] = 0;
        }

        return new MutableLongPoly[]{quotient.modulus(modulus), remainder.modulus(modulus)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDividerModulus(MutableLongPoly dividend, MutableLongPoly divider, long modulus, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = mod(-divider.cc(), modulus);
        long lcInverse = modInverse(divider.lc(), modulus);

        if (divider.lc() != 1)
            cc = mod(multiply(cc, lcInverse), modulus);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = mod(multiply(res, lcInverse), modulus);
            res = addMod(multiply(res, cc), tmp, modulus);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDivider(MutableLongPoly dividend, MutableLongPoly divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] pseudoDivideAndRemainderLinearDivider(MutableLongPoly dividend, MutableLongPoly divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, pow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] divideAndRemainderLinearDivider0(MutableLongPoly dividend, MutableLongPoly divider, long raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc();
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(raiseFactor, tmp));
            if (i == 0) break;
            if (res % lc != 0) return null;
            res = res / lc;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutableLongPoly[] pseudoDivideAndRemainderLinearDividerAdaptive(MutableLongPoly dividend, MutableLongPoly divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc(), factor = 1;
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(factor, tmp));
            if (i == 0) break;
            if (res % lc != 0) {
                long gcd = gcd(res, lc), f = lc / gcd;
                factor = multiply(factor, f);
                res = multiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = multiply(quotient[j], f);
            }
            res = res / lc;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutableLongPoly[]{MutableLongPoly.create(quotient), MutableLongPoly.create(res)};
    }


    /**
     * Returns the remainder of {@code dividend} divided by {@code divider} modulo {@code modulus}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  prime modulus
     * @return the remainder
     */
    public static MutableLongPoly remainder(final MutableLongPoly dividend,
                                            final MutableLongPoly divider,
                                            final long modulus,
                                            final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend).modulus(modulus);
        if (divider.degree == 0)
            return MutableLongPoly.zero();
        if (divider.degree == 1)
            return MutableLongPoly.create(dividend.evaluate(multiplyMod(-divider.cc(), modInverse(divider.lc(), modulus), modulus), modulus));
        return remainder0(dividend, divider, modulus, copy);
    }

    /** Plain school implementation */
    static MutableLongPoly remainder0(final MutableLongPoly dividend,
                                      final MutableLongPoly divider,
                                      final long modulus,
                                      final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly remainder = copy ? dividend.clone() : dividend;
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, divideMod(remainder.lc(), divider.lc(), modulus), i, modulus);

        return remainder.modulus(modulus);
    }

    /**
     * Returns the remainder of {@code dividend} divided by {@code divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return the remainder
     */
    public static MutableLongPoly remainder(final MutableLongPoly dividend,
                                            final MutableLongPoly divider,
                                            final boolean copy) {
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return MutableLongPoly.zero();
        if (divider.degree == 1) {
            if (divider.cc() % divider.lc() != 0)
                return null;
            return MutableLongPoly.create(dividend.evaluate(-divider.cc() / divider.lc()));
        }
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static MutableLongPoly remainder0(final MutableLongPoly dividend,
                                      final MutableLongPoly divider,
                                      final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutableLongPoly remainder = copy ? dividend.clone() : dividend;
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0)
                    return null;
                remainder.subtract(divider, remainder.lc() / divider.lc(), i);
            }
        return remainder;
    }
}
