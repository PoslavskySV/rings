package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class LinearAlgebra {
    private LinearAlgebra() {}

    public static BigInteger[] gaussianElimination(Domain<BigInteger> domain, BigInteger[][] A, BigInteger[] b) {
        assert A.length == b.length;
        int N = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (domain.compare(A[i][p], A[max][p]) > 0) {
                    max = i;
                }
            }
            BigInteger[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            BigInteger t = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (domain.isZero(A[p][p])) {
                throw new RuntimeException("Matrix is singular or nearly singular: " + Arrays.deepToString(A) + " = " + Arrays.toString(b));
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                BigInteger alpha = domain.divideAndRemainder(A[i][p], A[p][p])[0];
                b[i] = domain.subtract(b[i], domain.multiply(alpha, b[p]));
                for (int j = p; j < N; j++) {
                    A[i][j] = domain.subtract(A[i][j], domain.multiply(alpha, A[p][j]));
                }
            }
        }

        // back substitution
        BigInteger[] x = new BigInteger[N];
        for (int i = N - 1; i >= 0; i--) {
            BigInteger sum = domain.getZero();
            for (int j = i + 1; j < N; j++) {
                sum = domain.add(sum, domain.multiply(A[i][j], x[j]));
            }
            x[i] = domain.divideAndRemainder(domain.subtract(b[i], sum), A[i][i])[0];
        }
        return x;
    }

//    public static BigInteger vandermondeElimination(BigInteger[][] matrix, BigInteger[] rhs){
//
//    }
}
