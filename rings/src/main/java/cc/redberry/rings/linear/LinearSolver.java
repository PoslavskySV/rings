package cc.redberry.rings.linear;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.univar.UnivariateDivision;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

import static cc.redberry.rings.linear.LinearSolver.SystemInfo.*;

/**
 * Solver for quadratic linear system
 *
 * @since 1.0
 */
public final class LinearSolver {
    private LinearSolver() {}

    /**
     * Transpose square matrix
     */
    public static void transposeSquare(Object[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < i; ++j) {
                Object tmp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = tmp;
            }
        }
    }

    /**
     * Transpose square matrix
     */
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
     * Gives the row echelon form of the matrix
     *
     * @param ring   the ring
     * @param matrix the matrix
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Ring<E> ring, E[][] matrix) {
        return rowEchelonForm(ring, matrix, null, false, false);
    }

    /**
     * Gives the row echelon form of the matrix
     *
     * @param ring   the ring
     * @param matrix the matrix
     * @param reduce whether to calculate reduced row echelon form
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Ring<E> ring, E[][] matrix, boolean reduce) {
        return rowEchelonForm(ring, matrix, null, reduce, false);
    }

    /**
     * Gives the row echelon form of the linear system {@code lhs.x = rhs}.
     *
     * @param ring the ring
     * @param lhs  the lhs of the system
     * @param rhs  the rhs of the system
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Ring<E> ring, E[][] lhs, E[] rhs) {
        return rowEchelonForm(ring, lhs, rhs, false, false);
    }

    /**
     * Gives the row echelon form of the linear system {@code lhs.x = rhs}.
     *
     * @param ring                   the ring
     * @param lhs                    the lhs of the system
     * @param rhs                    the rhs of the system
     * @param reduce                 whether to calculate reduced row echelon form
     * @param breakOnUnderDetermined whether to return immediately if it was detected that system is under determined
     * @return the number of free variables
     */
    public static <E> int rowEchelonForm(Ring<E> ring, E[][] lhs, E[] rhs, boolean reduce, boolean breakOnUnderDetermined) {
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
            if (ring.isZero(lhs[row][iColumn])) {
                for (int iRow = row + 1; iRow < nRows; ++iRow)
                    if (!ring.isZero(lhs[iRow][iColumn])) {
                        max = iRow;
                        break;
                    }

                ArraysUtil.swap(lhs, row, max);
                if (rhs != null)
                    ArraysUtil.swap(rhs, row, max);
            }

            // singular
            if (ring.isZero(lhs[row][iColumn])) {
                if (breakOnUnderDetermined)
                    return 1;
                //nothing to do on this column
                ++nZeroColumns;
                to = Math.min(nRows + nZeroColumns, nColumns);
                continue;
            }

            // pivot within A and b
            for (int iRow = row + 1; iRow < nRows; ++iRow) {
                E alpha = ring.divideExact(lhs[iRow][iColumn], lhs[row][iColumn]);
                if (rhs != null)
                    rhs[iRow] = ring.subtract(rhs[iRow], ring.multiply(alpha, rhs[row]));
                if (!ring.isZero(alpha))
                    for (int iCol = iColumn; iCol < nColumns; ++iCol)
                        lhs[iRow][iCol] = ring.subtract(lhs[iRow][iCol], ring.multiply(alpha, lhs[row][iCol]));
            }
        }
        if (reduce)
            reducedRowEchelonForm(ring, lhs, rhs);
        return nZeroColumns;
    }

    /**
     * Gives the reduced row echelon form of the linear system {@code lhs.x = rhs} from a given row echelon form.
     *
     * @param ring the ring
     * @param lhs  the lhs of the system in the row echelon form
     * @param rhs  the rhs of the system
     */
    public static <E> void reducedRowEchelonForm(Ring<E> ring, E[][] lhs, E[] rhs) {
        int nRows = lhs.length;
        int nColumns = lhs[0].length;

        //number of zero columns
        int nZeroColumns = 0;
        for (int iColumn = 0, to = Math.min(nRows, nColumns); iColumn < to; ++iColumn) {
            // find pivot row and swap
            int iRow = iColumn - nZeroColumns;
            if (ring.isZero(lhs[iRow][iColumn])) {
                ++nZeroColumns;
                to = Math.min(nRows + nZeroColumns, nColumns);
                continue;
            }

            // scale current row
            E[] row = lhs[iRow];
            E val = row[iColumn];
            E valInv = ring.reciprocal(val);

            for (int i = iColumn; i < nColumns; i++)
                row[i] = ring.multiply(valInv, row[i]);
            if (rhs != null)
                rhs[iRow] = ring.multiply(valInv, rhs[iRow]);

            // scale all rows before
            for (int i = 0; i < iRow; i++) {
                E[] pRow = lhs[i];
                E v = pRow[iColumn];
                if (ring.isZero(v))
                    continue;
                for (int j = iColumn; j < nColumns; ++j)
                    pRow[j] = ring.subtract(pRow[j], ring.multiply(v, row[j]));
                if (rhs != null)
                    rhs[i] = ring.subtract(rhs[i], ring.multiply(v, rhs[iColumn]));
            }
        }
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and reduces lhs to row echelon form.
     *
     * @param ring the ring
     * @param lhs  the lhs of the system (will be reduced to row echelon form)
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solve(Ring<E> ring, E[][] lhs, E[] rhs) {
        int nUnknowns = lhs[0].length;
        if (nUnknowns == 0)
            return ring.createArray(0);
        E[] result = ring.createArray(nUnknowns);
        SystemInfo info = solve(ring, lhs, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Info about linear system
     */
    public enum SystemInfo {
        /** Under-determined system */
        UnderDetermined,
        /** Inconsistent system */
        Inconsistent,
        /** Consistent system */
        Consistent;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and reduces the lhs to row echelon form. The result is stored in {@code
     * result} (which should be of the enough length).
     *
     * @param ring   the ring
     * @param lhs    the lhs of the system (will be reduced to row echelon form)
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solve(Ring<E> ring, E[][] lhs, E[] rhs, E[] result) {
        if (lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            if (lhs[0].length == 1) {
                result[0] = ring.divideExact(rhs[0], lhs[0][0]);
                return Consistent;
            }
            if (lhs[0].length > 1)
                return UnderDetermined;
            if (lhs[0].length < 1)
                return Inconsistent;
        }

        int nUnderDetermined = rowEchelonForm(ring, lhs, rhs, false, true);
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
                if (!ring.isZero(rhs[i]))
                    // inconsistent system
                    return Inconsistent;


        // back substitution
        for (int i = nColumns - 1; i >= 0; i--) {
            E sum = ring.getZero();
            for (int j = i + 1; j < nColumns; j++) {
                sum = ring.add(sum, ring.multiply(lhs[i][j], result[j]));
            }
            result[i] = ring.divideExact(ring.subtract(rhs[i], sum), lhs[i][i]);
        }
        return Consistent;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result} (which should be of the enough
     * length).
     *
     * @param ring   the ring
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solve(Ring<E> ring, ArrayList<E[]> lhs, ArrayList<E> rhs, E[] result) {
        return solve(ring, lhs.toArray(ring.createArray2d(lhs.size())), rhs.toArray(ring.createArray(rhs.size())), result);
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 +
     * ... row[i]^N * xN = rhs[i] }).
     *
     * @param ring the ring
     * @param row  the Vandermonde coefficients
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solveVandermonde(Ring<E> ring, E[] row, E[] rhs) {
        E[] result = ring.createArray(rhs.length);
        SystemInfo info = solveVandermonde(ring, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +
     * row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }).
     *
     * @param ring the ring
     * @param row  the Vandermonde coefficients
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static <E> E[] solveVandermondeT(Ring<E> ring, E[] row, E[] rhs) {
        E[] result = ring.createArray(rhs.length);
        SystemInfo info = solveVandermondeT(ring, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 +
     * ... row[i]^N * xN = rhs[i] }) and stores the result in {@code result} (which should be of the enough length).
     *
     * @param ring   the ring
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solveVandermonde(Ring<E> ring, E[] row, E[] rhs, E[] result) {
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
        UnivariatePolynomial<E> master = UnivariatePolynomial.one(ring);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(ring.negate(row[i]), ring.getOne());
            master = master.multiply(lins[i]);
        }


        for (int i = 0; i < result.length; i++)
            result[i] = ring.getZero();

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomial<E> quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            E cf = quot.evaluate(row[i]);
            if (ring.isZero(cf))
                return UnderDetermined;
            quot = quot.divideOrNull(cf);
            if (quot == null)
                throw new IllegalArgumentException();
            for (int j = 0; j < row.length; ++j)
                result[j] = ring.add(result[j], ring.multiply(rhs[i], quot.get(j)));
        }
        return Consistent;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +
     * row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }) and stores the result in {@code result} (which should be of the
     * enough length).
     *
     * @param ring   the ring
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static <E> SystemInfo solveVandermondeT(Ring<E> ring, E[] row, E[] rhs, E[] result) {
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
        UnivariatePolynomial<E> master = UnivariatePolynomial.one(ring);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(ring.negate(row[i]), ring.getOne());
            master = master.multiply(lins[i]);
        }

        for (int i = 0; i < row.length; i++) {
            UnivariatePolynomial<E> quot = UnivariateDivision.divideAndRemainder(master, lins[i], true)[0];
            E cf = quot.evaluate(row[i]);
            if (ring.isZero(cf))
                return UnderDetermined;
            quot = quot.divideOrNull(cf);
            if (quot == null)
                throw new IllegalArgumentException();
            result[i] = ring.getZero();
            for (int j = 0; j < row.length; ++j)
                result[i] = ring.add(result[i], ring.multiply(rhs[j], quot.get(j)));
        }
        return Consistent;
    }


    /* ========================================= Machine numbers ============================================ */

    /**
     * Gives the row echelon form of the matrix
     *
     * @param ring   the ring
     * @param matrix the matrix
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 ring, long[][] matrix) {
        return rowEchelonForm(ring, matrix, false);
    }

    /**
     * Gives the row echelon form of the matrix
     *
     * @param ring   the ring
     * @param matrix the matrix
     * @param reduce whether to calculate reduced row echelon form
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 ring, long[][] matrix, boolean reduce) {
        return rowEchelonForm(ring, matrix, null, reduce, false);
    }

    /**
     * Gives the row echelon form of the linear system {@code lhs.x = rhs} (rhs may be null).
     *
     * @param ring the ring
     * @param lhs  the lhs of the system
     * @param rhs  the rhs of the system (may be null)
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 ring, long[][] lhs, long[] rhs) {
        return rowEchelonForm(ring, lhs, rhs, false, false);
    }

    /**
     * Gives the row echelon form of the linear system {@code lhs.x = rhs} (rhs may be null).
     *
     * @param ring                   the ring
     * @param lhs                    the lhs of the system
     * @param rhs                    the rhs of the system (may be null)
     * @param reduce                 whether to calculate reduced row echelon form
     * @param breakOnUnderDetermined whether to return immediately if it was detected that system is under determined
     * @return the number of free variables
     */
    public static int rowEchelonForm(IntegersZp64 ring, long[][] lhs, long[] rhs, boolean reduce, boolean breakOnUnderDetermined) {
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
            int nonZero = row;
            if (lhs[row][iColumn] == 0) {
                for (int iRow = row + 1; iRow < nRows; ++iRow)
                    if (lhs[iRow][iColumn] != 0) {
                        nonZero = iRow;
                        break;
                    }

                ArraysUtil.swap(lhs, row, nonZero);
                if (rhs != null)
                    ArraysUtil.swap(rhs, row, nonZero);
            }

            // singular
            if (lhs[row][iColumn] == 0) {
                if (breakOnUnderDetermined)
                    return 1;
                //nothing to do on this column
                ++nZeroColumns;
                to = Math.min(nRows + nZeroColumns, nColumns);
                continue;
            }

            // pivot within A and b
            for (int iRow = row + 1; iRow < nRows; ++iRow) {
                long alpha = ring.divide(lhs[iRow][iColumn], lhs[row][iColumn]);
                if (rhs != null)
                    rhs[iRow] = ring.subtract(rhs[iRow], ring.multiply(alpha, rhs[row]));
                if (alpha != 0)
                    for (int iCol = iColumn; iCol < nColumns; ++iCol)
                        lhs[iRow][iCol] = ring.subtract(lhs[iRow][iCol], ring.multiply(alpha, lhs[row][iCol]));
            }
        }
        if (reduce)
            reducedRowEchelonForm(ring, lhs, rhs);
        return nZeroColumns;
    }

    /**
     * Gives the reduced row echelon form of the linear system {@code lhs.x = rhs} from a given row echelon form.
     *
     * @param ring the ring
     * @param lhs  the lhs of the system in the row echelon form
     * @param rhs  the rhs of the system
     */
    public static void reducedRowEchelonForm(IntegersZp64 ring, long[][] lhs, long[] rhs) {
        int nRows = lhs.length;
        int nColumns = lhs[0].length;

        //number of zero columns
        int nZeroColumns = 0;
        for (int iColumn = 0, to = Math.min(nRows, nColumns); iColumn < to; ++iColumn) {
            // find pivot row and swap
            int iRow = iColumn - nZeroColumns;
            if (lhs[iRow][iColumn] == 0) {
                ++nZeroColumns;
                to = Math.min(nRows + nZeroColumns, nColumns);
                continue;
            }

            // scale current row
            long[] row = lhs[iRow];
            long val = row[iColumn];
            long valInv = ring.reciprocal(val);

            for (int i = iColumn; i < nColumns; i++)
                row[i] = ring.multiply(valInv, row[i]);
            if (rhs != null)
                rhs[iRow] = ring.multiply(valInv, rhs[iRow]);

            // scale all rows before
            for (int i = 0; i < iRow; i++) {
                long[] pRow = lhs[i];
                long v = pRow[iColumn];
                if (v == 0)
                    continue;
                for (int j = iColumn; j < nColumns; ++j)
                    pRow[j] = ring.subtract(pRow[j], ring.multiply(v, row[j]));
                if (rhs != null)
                    rhs[i] = ring.subtract(rhs[i], ring.multiply(v, rhs[iColumn]));
            }
        }
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and reduces the lhs to row echelon form.
     *
     * @param ring the ring
     * @param lhs  the lhs of the system  (will be reduced to row echelon form)
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solve(IntegersZp64 ring, long[][] lhs, long[] rhs) {
        int nUnknowns = lhs[0].length;
        if (nUnknowns == 0)
            return new long[0];
        long[] result = new long[nUnknowns];
        SystemInfo info = solve(ring, lhs, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and reduces the lhs to row echelon form. The result is stored in {@code
     * result} (which should be of the enough length).
     *
     * @param ring   the ring
     * @param lhs    the lhs of the system  (will be reduced to row echelon form)
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solve(IntegersZp64 ring, long[][] lhs, long[] rhs, long[] result) {
        if (lhs.length != rhs.length)
            throw new IllegalArgumentException("lhs.length != rhs.length");
        if (rhs.length == 0)
            return Consistent;
        if (rhs.length == 1) {
            if (lhs[0].length == 1) {
                result[0] = ring.divide(rhs[0], lhs[0][0]);
                return Consistent;
            }
            if (lhs[0].length > 1)
                return UnderDetermined;
            if (lhs[0].length < 1)
                return Inconsistent;
        }

        int nUnderDetermined = rowEchelonForm(ring, lhs, rhs, false, true);
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
                sum = ring.add(sum, ring.multiply(lhs[i][j], result[j]));
            }
            result[i] = ring.divide(ring.subtract(rhs[i], sum), lhs[i][i]);
        }
        return Consistent;
    }

    /**
     * Solves linear system {@code lhs.x = rhs} and stores the result in {@code result} (which should be of the enough
     * length).
     *
     * @param ring   the ring
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solve(IntegersZp64 ring, ArrayList<long[]> lhs, TLongArrayList rhs, long[] result) {
        return solve(ring, lhs.toArray(new long[lhs.size()][]), rhs.toArray(), result);
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 +
     * ... row[i]^N * xN = rhs[i] }).
     *
     * @param ring the ring
     * @param row  the Vandermonde coefficients
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solveVandermonde(IntegersZp64 ring, long[] row, long[] rhs) {
        long[] result = new long[rhs.length];
        SystemInfo info = solveVandermonde(ring, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +
     * row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }).
     *
     * @param ring the ring
     * @param row  the Vandermonde coefficients
     * @param rhs  the rhs of the system
     * @return the solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static long[] solveVandermondeT(IntegersZp64 ring, long[] row, long[] rhs) {
        long[] result = new long[rhs.length];
        SystemInfo info = solveVandermondeT(ring, row, rhs, result);
        if (info != Consistent)
            throw new ArithmeticException("singular or under-determined matrix");
        return result;
    }

    /**
     * Solves Vandermonde linear system (that is with i-th equation of the form {@code row[i]^0 * x0 +  row[i]^1 * x1 +
     * ... row[i]^N * xN = rhs[i] }) and stores the result in {@code result} (which should be of the enough length).
     *
     * @param ring   the ring
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solveVandermonde(IntegersZp64 ring, long[] row, long[] rhs, long[] result) {
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
        UnivariatePolynomialZp64 master = UnivariatePolynomialZp64.one(ring);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(ring.negate(row[i]), 1L);
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
                result[j] = ring.add(result[j], ring.multiply(rhs[i], quot.get(j)));
        }
        return Consistent;
    }

    /**
     * Solves transposed Vandermonde linear system (that is with i-th equation of the form {@code row[0]^i * x0 +
     * row[1]^i * x1 + ... row[N]^i * xN = rhs[i] }) and stores the result in {@code result} (which should be of the
     * enough length).
     *
     * @param ring   the ring
     * @param row    the Vandermonde coefficients
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solveVandermondeT(IntegersZp64 ring, long[] row, long[] rhs, long[] result) {
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
        UnivariatePolynomialZp64 master = UnivariatePolynomialZp64.one(ring);
        for (int i = 0; i < row.length; ++i) {
            lins[i] = master.createLinear(ring.negate(row[i]), 1L);
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
                result[i] = ring.add(result[i], ring.multiply(rhs[j], quot.get(j)));
        }
        return Consistent;
    }
}
