package cc.r2.core.poly.multivar;

import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersZp64;
import cc.r2.core.poly.univar.UnivariateDivision;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariatePolynomialZp64;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

import static cc.r2.core.poly.multivar.LinearAlgebra.SystemInfo.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class LinearAlgebra {
    private LinearAlgebra() {}

    public static void transposeSquare(Object[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < i; ++j) {
                Object tmp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = tmp;
            }
        }
    }

    public static void transposeSquare(long[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < i; ++j) {
                long tmp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = tmp;
            }
        }
    }

    /**
     * Gives the row‐reduced form of the matrix
     *
     * @param domain the domain
     * @param matrix the matrix
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Domain<E> domain, E[][] matrix) {
        return rowEchelonForm(domain, matrix, null);
    }

    /**
     * Gives the row‐reduced form of the linear system {@code lhs.x = rhs}.
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Domain<E> domain, E[][] lhs, E[] rhs) {
        if (rhs != null && lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");

        if (lhs.length == 0)
            return 0;

        int nRows = lhs.length;
        int nColumns = lhs[0].length;

        //number of zero columns
        int nZeroColumns = 0;
        for (int iColumn = 0, to = Math.min(nRows, nColumns); iColumn < to; ++iColumn) {

            // find pivot row and swap
            int row = iColumn - nZeroColumns;
            int max = row;
            for (int iRow = row + 1; iRow < nRows; ++iRow)
                if (domain.compare(lhs[iRow][iColumn], lhs[max][iColumn]) > 0)
                    max = iRow;

            ArraysUtil.swap(lhs, row, max);
            if (rhs != null)
                ArraysUtil.swap(rhs, row, max);

            // singular
            if (domain.isZero(lhs[row][iColumn])) {
                //nothing to do on this column
                ++nZeroColumns;
                continue;
            }

            // pivot within A and b
            for (int iRow = row + 1; iRow < nRows; ++iRow) {
                E alpha = domain.divideExact(lhs[iRow][iColumn], lhs[row][iColumn]);
                if (rhs != null)
                    rhs[iRow] = domain.subtract(rhs[iRow], domain.multiply(alpha, rhs[row]));
                if (!domain.isZero(alpha))
                    for (int iCol = iColumn; iCol < nColumns; ++iCol)
                        lhs[iRow][iCol] = domain.subtract(lhs[iRow][iCol], domain.multiply(alpha, lhs[row][iCol]));
            }
        }
        return nZeroColumns;
    }

    /**
     * Solves linear system {@code lhs.x = rhs}.
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solve(Domain<E> domain, E[][] lhs, E[] rhs) {
        int nUnknowns = lhs[0].length;
        if (nUnknowns == 0)
            return domain.createArray(0);
        E[] result = domain.createArray(nUnknowns);
        SystemInfo info = solve(domain, lhs, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    enum SystemInfo {
        /** Under-determined system */
        UnderDetermined,
        /** Inconsistent system */
        Inconsistent,
        /** Consistent system */
        Consistent
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solve(Domain<E> domain, E[][] lhs, E[] rhs, E[] result) {
        if (lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            if (lhs[0].length == 1) {
                result[0] = domain.divideExact(rhs[0], lhs[0][0]);
                return Consistent;
            }
            if (lhs[0].length > 1)
                return UnderDetermined;
            if (lhs[0].length < 1)
                return Inconsistent;
        }

        int nUnderDetermined = rowEchelonForm(domain, lhs, rhs);
        if (nUnderDetermined > 0)
            // under-determined system
            return UnderDetermined;

        int nRows = rhs.length;
        int nColumns = lhs[0].length;

        if (nColumns > nRows)
            // under-determined system
            return null;

        if (nRows > nColumns)
            // over-determined system
            // check that all rhs are zero
            for (int i = nColumns; i < nRows; ++i)
                if (!domain.isZero(rhs[i]))
                    // inconsistent system
                    return Inconsistent;


        // back substitution
        for (int i = nColumns - 1; i >= 0; i--) {
            E sum = domain.getZero();
            for (int j = i + 1; j < nColumns; j++) {
                sum = domain.add(sum, domain.multiply(lhs[i][j], result[j]));
            }
            result[i] = domain.divideExact(domain.subtract(rhs[i], sum), lhs[i][i]);
        }
        return Consistent;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solve(Domain<E> domain, ArrayList<E[]> lhs, ArrayList<E> rhs, E[] result) {
        return solve(domain, lhs.toArray(domain.createArray2d(lhs.size())), rhs.toArray(domain.createArray(rhs.size())), result);
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 + ... row[i]^N * xN = rhs[i] }).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solveVandermonde(Domain<E> domain, E[] row, E[] rhs) {
        E[] result = domain.createArray(rhs.length);
        SystemInfo info = solveVandermonde(domain, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +  row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solveVandermondeT(Domain<E> domain, E[] row, E[] rhs) {
        E[] result = domain.createArray(rhs.length);
        SystemInfo info = solveVandermondeT(domain, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 + ... row[i]^N * xN = rhs[i] })
     * and stores the result in {@code result} (which should be of the enough length).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solveVandermonde(Domain<E> domain, E[] row, E[] rhs, E[] result) {
        if (row.length != rhs.length)
            throw new IllegalArgumentException("not a square Vandermonde matrix");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            result[0] = rhs[0];
            return Consistent;
        }
        @SuppressWarnings("unchecked")
        UnivariatePolynomial<E>[] lins = new UnivariatePolynomial[row.length];
        UnivariatePolynomial<E> master = UnivariatePolynomial.one(domain);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(domain.negate(row[i]), domain.getOne());
            master = master.multiply(lins[i]);
        }


        for (int i = 0; i < result.length; i++)
            result[i] = domain.getZero();

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomial<E> quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            E cf = quot.evaluate(row[i]);
            if (domain.isZero(cf))
                return UnderDetermined;
            quot = quot.divideOrNull(cf);
            if (quot == null)
                throw new IllegalArgumentException();
            for (int j = 0; j < row.length; ++j)
                result[j] = domain.add(result[j], domain.multiply(rhs[i], quot.get(j)));
        }
        return Consistent;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +  row[1]^i * x1 + ... row[N]^i * xN = rhs[i] })
     * and stores the result in {@code result} (which should be of the enough length).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solveVandermondeT(Domain<E> domain, E[] row, E[] rhs, E[] result) {
        if (row.length != rhs.length)
            throw new IllegalArgumentException("not a square Vandermonde matrix");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            result[0] = rhs[0];
            return Consistent;
        }
        @SuppressWarnings("unchecked")
        UnivariatePolynomial<E>[] lins = new UnivariatePolynomial[row.length];
        UnivariatePolynomial<E> master = UnivariatePolynomial.one(domain);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(domain.negate(row[i]), domain.getOne());
            master = master.multiply(lins[i]);
        }

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomial<E> quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            E cf = quot.evaluate(row[i]);
            if (domain.isZero(cf))
                return UnderDetermined;
            quot = quot.divideOrNull(cf);
            if (quot == null)
                throw new IllegalArgumentException();
            result[i] = domain.getZero();
            for (int j = 0; j < row.length; ++j)
                result[i] = domain.add(result[i], domain.multiply(rhs[j], quot.get(j)));
        }
        return Consistent;
    }


    /* ========================================= Machine numbers ============================================ */

    /**
     * Gives the row‐reduced form of the matrix
     *
     * @param domain the domain
     * @param matrix the matrix
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 domain, long[][] matrix) {
        return rowEchelonForm(domain, matrix, null);
    }

    /**
     * Gives the row‐reduced form of the linear system {@code lhs.x = rhs}.
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 domain, long[][] lhs, long[] rhs) {
        if (rhs != null && lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");

        if (lhs.length == 0)
            return 0;

        int nRows = lhs.length;
        int nColumns = lhs[0].length;

        //number of zero columns
        int nZeroColumns = 0;
        for (int iColumn = 0, to = Math.min(nRows, nColumns); iColumn < to; ++iColumn) {

            // find pivot row and swap
            int row = iColumn - nZeroColumns;
            int max = row;
            for (int iRow = row + 1; iRow < nRows; ++iRow)
                if (lhs[iRow][iColumn] > lhs[max][iColumn])
                    max = iRow;

            ArraysUtil.swap(lhs, row, max);
            if (rhs != null)
                ArraysUtil.swap(rhs, row, max);

            // singular
            if (lhs[row][iColumn] == 0) {
                //nothing to do on this column
                ++nZeroColumns;
                continue;
            }

            // pivot within A and b
            for (int iRow = row + 1; iRow < nRows; ++iRow) {
                long alpha = domain.divide(lhs[iRow][iColumn], lhs[row][iColumn]);
                if (rhs != null)
                    rhs[iRow] = domain.subtract(rhs[iRow], domain.multiply(alpha, rhs[row]));
                if (alpha != 0)
                    for (int iCol = iColumn; iCol < nColumns; ++iCol)
                        lhs[iRow][iCol] = domain.subtract(lhs[iRow][iCol], domain.multiply(alpha, lhs[row][iCol]));
            }
        }
        return nZeroColumns;
    }

    /**
     * Solves linear system {@code lhs.x = rhs}.
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solve(IntegersZp64 domain, long[][] lhs, long[] rhs) {
        int nUnknowns = lhs[0].length;
        if (nUnknowns == 0)
            return new long[0];
        long[] result = new long[nUnknowns];
        SystemInfo info = solve(domain, lhs, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solve(IntegersZp64 domain, long[][] lhs, long[] rhs, long[] result) {
        if (lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            if (lhs[0].length == 1) {
                result[0] = domain.divide(rhs[0], lhs[0][0]);
                return Consistent;
            }
            if (lhs[0].length > 1)
                return UnderDetermined;
            if (lhs[0].length < 1)
                return Inconsistent;
        }

        int nUnderDetermined = rowEchelonForm(domain, lhs, rhs);
        if (nUnderDetermined > 0)
            // under-determined system
            return UnderDetermined;

        int nRows = rhs.length;
        int nColumns = lhs[0].length;

        if (nColumns > nRows)
            // under-determined system
            return null;

        if (nRows > nColumns)
            // over-determined system
            // check that all rhs are zero
            for (int i = nColumns; i < nRows; ++i)
                if (rhs[i] != 0)
                    // inconsistent system
                    return Inconsistent;


        // back substitution
        for (int i = nColumns - 1; i >= 0; i--) {
            long sum = 0;
            for (int j = i + 1; j < nColumns; j++) {
                sum = domain.add(sum, domain.multiply(lhs[i][j], result[j]));
            }
            result[i] = domain.divide(domain.subtract(rhs[i], sum), lhs[i][i]);
        }
        return Consistent;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solve(IntegersZp64 domain, ArrayList<long[]> lhs, TLongArrayList rhs, long[] result) {
        return solve(domain, lhs.toArray(new long[lhs.size()][]), rhs.toArray(), result);
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 + ... row[i]^N * xN = rhs[i] }).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solveVandermonde(IntegersZp64 domain, long[] row, long[] rhs) {
        long[] result = new long[rhs.length];
        SystemInfo info = solveVandermonde(domain, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +  row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solveVandermondeT(IntegersZp64 domain, long[] row, long[] rhs) {
        long[] result = new long[rhs.length];
        SystemInfo info = solveVandermondeT(domain, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 + ... row[i]^N * xN = rhs[i] })
     * and stores the result in {@code result} (which should be of the enough length).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solveVandermonde(IntegersZp64 domain, long[] row, long[] rhs, long[] result) {
        if (row.length != rhs.length)
            throw new IllegalArgumentException("not a square Vandermonde matrix");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            result[0] = rhs[0];
            return Consistent;
        }
        @SuppressWarnings("unchecked")
        UnivariatePolynomialZp64[] lins = new UnivariatePolynomialZp64[row.length];
        UnivariatePolynomialZp64 master = UnivariatePolynomialZp64.one(domain);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(domain.negate(row[i]), 1L);
            master = master.multiply(lins[i]);
        }


        for (int i = 0; i < result.length; i++)
            result[i] = 0;

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomialZp64 quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            long cf = quot.evaluate(row[i]);
            if (cf == 0)
                return UnderDetermined;
            quot = quot.divide(cf);
            for (int j = 0; j < row.length; ++j)
                result[j] = domain.add(result[j], domain.multiply(rhs[i], quot.get(j)));
        }
        return Consistent;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +  row[1]^i * x1 + ... row[N]^i * xN = rhs[i] })
     * and stores the result in {@code result} (which should be of the enough length).
     *
     * @param domain the domain
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solveVandermondeT(IntegersZp64 domain, long[] row, long[] rhs, long[] result) {
        if (row.length != rhs.length)
            throw new IllegalArgumentException("not a square Vandermonde matrix");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            result[0] = rhs[0];
            return Consistent;
        }
        @SuppressWarnings("unchecked")
        UnivariatePolynomialZp64[] lins = new UnivariatePolynomialZp64[row.length];
        UnivariatePolynomialZp64 master = UnivariatePolynomialZp64.one(domain);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(domain.negate(row[i]), 1L);
            master = master.multiply(lins[i]);
        }

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomialZp64 quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            long cf = quot.evaluate(row[i]);
            if (cf == 0)
                return UnderDetermined;
            quot = quot.divide(cf);
            if (quot == null)
                throw new IllegalArgumentException();
            result[i] = 0;
            for (int j = 0; j < row.length; ++j)
                result[i] = domain.add(result[i], domain.multiply(rhs[j], quot.get(j)));
        }
        return Consistent;
    }
}
