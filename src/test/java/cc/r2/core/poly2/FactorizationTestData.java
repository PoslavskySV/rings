package cc.r2.core.poly2;

import gnu.trove.list.array.TIntArrayList;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

public final class FactorizationTestData {
    public final MutablePolynomialMod poly;
    public final FactorDecomposition factorization;
    public final long modulus;

    public FactorizationTestData(MutablePolynomialMod poly, FactorDecomposition factorization, long modulus) {
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
        MutablePolynomialMod poly = MutablePolynomialZ.create(parseArray(parts[1])).modulus(modulus);
        MutablePolynomialMod[] factors = new MutablePolynomialMod[parts.length - 2];
        int[] exponents = new int[parts.length - 2];
        for (int i = 2; i < parts.length; i++) {
            long[] data = parseArray(parts[i].trim());
            exponents[i - 2] = (int) data[0];
            factors[i - 2] = MutablePolynomialZ.create(Arrays.copyOfRange(data, 1, data.length)).modulus(modulus);
        }
        return new FactorizationTestData(poly, new lFactorDecomposition(new ArrayList<>(Arrays.asList(factors)), new TIntArrayList(exponents), 1), modulus);
    }

    public static Object[] decodePolynomial(String string) {
        String[] coefficients = string.split(",");
        int exponent = Integer.parseInt(coefficients[0]);
        long[] data = new long[coefficients.length - 1];
        for (int i = 1; i < coefficients.length; i++)
            data[i - 1] = Long.parseLong(coefficients[i]);
        return new Object[]{MutablePolynomialZ.create(data), exponent};
    }

    public static FactorizationTestData decodeModFactorization(String string) {
        String[] parts = string.split("\\|");
        long modulus = Long.parseLong(parts[0]);
        Object[] poly = decodePolynomial(parts[1]);
        int[] exponents = new int[parts.length - 2];
        MutablePolynomialMod[] factors = new MutablePolynomialMod[parts.length - 2];
        for (int i = 2; i < parts.length; i++) {
            Object[] pe = decodePolynomial(parts[i]);
            factors[i - 2] = ((MutablePolynomialZ) pe[0]).modulus(modulus);
            exponents[i - 2] = (int) pe[1];
        }
        return new FactorizationTestData(
                PolynomialArithmetics.polyPow(((MutablePolynomialZ) poly[0]).modulus(modulus), (int) poly[1], false),
                new lFactorDecomposition(new ArrayList<>(Arrays.asList(factors)), new TIntArrayList(exponents), 1),
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
