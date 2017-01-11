package cc.r2.core.polynomial;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Iterator;

public final class FactorizationTestData {
    public final MutableLongPoly poly;
    public final SmallPolynomials.Factorization factorization;
    public final long modulus;

    public FactorizationTestData(MutableLongPoly poly, SmallPolynomials.Factorization factorization, long modulus) {
        this.poly = poly;
        this.factorization = factorization;
        this.modulus = modulus;
    }

    private static long[] parseArray(String string) {
        String[] split = string.split(",");
        long[] data = new long[split.length];
        for (int i = 0; i < split.length; i++)
            data[i] = Long.parseLong(split[i].trim());
        return data;
    }

    public static FactorizationTestData decode(String line) {
        String[] parts = line.split("\\|");
        long modulus = Long.parseLong(parts[0]);
        MutableLongPoly poly = MutableLongPoly.create(parseArray(parts[1]));
        MutableLongPoly[] factors = new MutableLongPoly[parts.length - 2];
        int[] exponents = new int[parts.length - 2];
        for (int i = 2; i < parts.length; i++) {
            long[] data = parseArray(parts[i].trim());
            exponents[i - 2] = (int) data[0];
            factors[i - 2] = MutableLongPoly.create(Arrays.copyOfRange(data, 1, data.length));
        }
        return new FactorizationTestData(poly, new SmallPolynomials.Factorization(factors, exponents, 1), modulus);
    }

    public static Object[] decodePolynomial(String string) {
        String[] coefficients = string.split(",");
        int exponent = Integer.parseInt(coefficients[0]);
        long[] data = new long[coefficients.length - 1];
        for (int i = 1; i < coefficients.length; i++)
            data[i - 1] = Long.parseLong(coefficients[i]);
        return new Object[]{MutableLongPoly.create(data), exponent};
    }

    public static FactorizationTestData decodeModFactorization(String string) {
        String[] parts = string.split("\\|");
        long modulus = Long.parseLong(parts[0]);
        Object[] poly = decodePolynomial(parts[1]);
        int[] exponents = new int[parts.length - 2];
        MutableLongPoly[] factors = new MutableLongPoly[parts.length - 2];
        for (int i = 2; i < parts.length; i++) {
            Object[] pe = decodePolynomial(parts[i]);
            factors[i - 2] = (MutableLongPoly) pe[0];
            exponents[i - 2] = (int) pe[1];
        }
        return new FactorizationTestData(
                SmallPolynomialArithmetics.polyPowMod((MutableLongPoly) poly[0], (int) poly[1], modulus, false),
                new SmallPolynomials.Factorization(factors, exponents, 1),
                modulus);
    }

    public static Iterable<FactorizationTestData> allMod(InputStream is) {
        return () -> allModIterator(is);
    }

    public static Iterator<FactorizationTestData> allModIterator(InputStream is) {
        final BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        return new Iterator<FactorizationTestData>() {
            String next = null;

            @Override
            public boolean hasNext() {
                try {
                    return (next = reader.readLine()) != null;
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }

            @Override
            public FactorizationTestData next() {
                return next == null ? null : decodeModFactorization(next);
            }
        };
    }
}
