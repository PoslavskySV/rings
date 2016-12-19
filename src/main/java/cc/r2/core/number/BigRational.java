/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package cc.r2.core.number;

import java.io.Serializable;

public final class BigRational
        implements FieldElement<BigRational>, Serializable {
    /**
     * A fraction representing "2 / 1".
     */
    public static final BigRational TWO = new BigRational(2);

    /**
     * A fraction representing "1".
     */
    public static final BigRational ONE = new BigRational(1);

    /**
     * A fraction representing "0".
     */
    public static final BigRational ZERO = new BigRational(0);

    /**
     * A fraction representing "-1 / 1".
     */
    public static final BigRational MINUS_ONE = new BigRational(-1);

    /**
     * A fraction representing "4/5".
     */
    public static final BigRational FOUR_FIFTHS = new BigRational(4, 5);

    /**
     * A fraction representing "1/5".
     */
    public static final BigRational ONE_FIFTH = new BigRational(1, 5);

    /**
     * A fraction representing "1/2".
     */
    public static final BigRational ONE_HALF = new BigRational(1, 2);

    /**
     * A fraction representing "1/4".
     */
    public static final BigRational ONE_QUARTER = new BigRational(1, 4);

    /**
     * A fraction representing "1/3".
     */
    public static final BigRational ONE_THIRD = new BigRational(1, 3);

    /**
     * A fraction representing "3/5".
     */
    public static final BigRational THREE_FIFTHS = new BigRational(3, 5);

    /**
     * A fraction representing "3/4".
     */
    public static final BigRational THREE_QUARTERS = new BigRational(3, 4);

    /**
     * A fraction representing "2/5".
     */
    public static final BigRational TWO_FIFTHS = new BigRational(2, 5);

    /**
     * A fraction representing "2/4".
     */
    public static final BigRational TWO_QUARTERS = new BigRational(2, 4);

    /**
     * A fraction representing "2/3".
     */
    public static final BigRational TWO_THIRDS = new BigRational(2, 3);

    /**
     * Serializable version identifier.
     */
    private static final long serialVersionUID = -5630213147331578515L;

    /**
     * <code>BigInteger</code> representation of 100.
     */
    private static final BigInteger ONE_HUNDRED = BigInteger.valueOf(100);

    /**
     * The numerator.
     */
    private final BigInteger numerator;

    /**
     * The denominator.
     */
    private final BigInteger denominator;

    public BigRational(int num) {
        this(BigInteger.valueOf(num), BigInteger.ONE);
    }

    public BigRational(int num, int den) {
        this(BigInteger.valueOf(num), BigInteger.valueOf(den));
    }

    public BigRational(BigInteger num, BigInteger den) {
        if (den.isZero())
            throw new IllegalArgumentException("Denominator is zero");
        if (num.isZero()) {
            numerator = BigInteger.ZERO;
            denominator = BigInteger.ONE;
        } else {

            // reduce numerator and denominator by greatest common denominator
            final BigInteger gcd = num.gcd(den);
            if (BigInteger.ONE.compareTo(gcd) < 0) {
                num = num.divide(gcd);
                den = den.divide(gcd);
            }

            // move sign to numerator
            if (den.signum() < 0) {
                num = num.negate();
                den = den.negate();
            }

            // store the values in the final fields
            numerator = num;
            denominator = den;

        }
    }

    private BigRational(boolean nogcd, BigInteger num, BigInteger den) {
        // move sign to numerator
        if (den.signum() < 0) {
            num = num.negate();
            den = den.negate();
        }

        this.numerator = num;
        this.denominator = den;
    }

    public BigInteger getDenominator() {
        return denominator;
    }

    public BigInteger getNumerator() {
        return numerator;
    }

    @Override
    public Ring<BigRational> getRing() {
        return BigRationalField.BigRationalField;
    }

    public BigRational abs() {
        return (numerator.signum() >= 0) ? this : negate();
    }

    @Override
    public BigRational add(final BigRational fraction) {
        if (fraction.isZero())
            return this;

        BigInteger num, den;
        if (denominator.equals(fraction.denominator)) {
            num = numerator.add(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).add((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return new BigRational(num, den);

    }

    @Override
    public BigRational multiply(final BigRational fraction) {
        if (numerator.isZero() || fraction.numerator.isZero())
            return ZERO;
        return new BigRational(numerator.multiply(fraction.numerator), denominator.multiply(fraction.denominator));
    }

    @Override
    public BigRational divide(final BigRational fraction) {
        if (fraction.numerator.isZero())
            throw new IllegalArgumentException("Division by zero");

        return multiply(fraction.reciprocal());
    }

    @Override
    public BigRational negate() {
        if (isZero())
            return this;
        return new BigRational(true, numerator.negate(), denominator);
    }

    @Override
    public BigRational reciprocal() {
        return new BigRational(true, denominator, numerator);
    }

    public BigRational subtract(final BigRational fraction) {
        if (fraction.isZero())
            return this;

        BigInteger num, den;
        if (denominator.equals(fraction.denominator)) {
            num = numerator.subtract(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).subtract((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return new BigRational(num, den);
    }

    @Override
    public boolean isZero() {
        return numerator.isZero();
    }

    @Override
    public boolean isOne() {
        return numerator.isOne() && denominator.isOne();
    }

    @Override
    public BigRational gcd(BigRational oth) {
        return new BigRational(numerator.multiply(oth.denominator).gcd(oth.numerator.multiply(denominator)), denominator.multiply(oth.denominator));
    }

    @Override
    public BigRational[] divideAndRemainder(BigRational oth) {
        return new BigRational[]{this.divide(oth), ZERO};
    }

    @Override
    public BigRational getZero() {
        return ZERO;
    }

    @Override
    public BigRational getOne() {
        return ONE;
    }

    @Override
    public int compareTo(final BigRational object) {
        BigInteger nOd = numerator.multiply(object.denominator);
        BigInteger dOn = denominator.multiply(object.numerator);
        return nOd.compareTo(dOn);
    }


    @Override
    public boolean equals(final Object other) {
        if (this == other)
            return true;

        if (!(other instanceof BigRational))
            return false;
        BigRational rhs = (BigRational) other;
        return numerator.equals(rhs.numerator) && denominator.equals(rhs.denominator);
    }

    @Override
    public int hashCode() {
        return 37 * (37 * 17 + numerator.hashCode()) + denominator.hashCode();
    }

    @Override
    public String toString() {
        String str;
        if (denominator.isOne()) {
            str = numerator.toString();
        } else if (numerator.isZero()) {
            str = "0";
        } else {
            str = numerator + "/" + denominator;
        }
        return str;
    }
}
