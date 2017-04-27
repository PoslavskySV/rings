package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.util.ArraysUtil;

import java.util.ArrayList;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class LinearAlgebra {
    private LinearAlgebra() {}

    /**
     * Gives the row‐reduced form of the matrix
     *
     * @param domain the domain
     * @param matrix the matrix
     * @return the number of free variables
     */
    public static int rowEchelonForm(Domain<BigInteger> domain, BigInteger[][] matrix) {
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
    public static int rowEchelonForm(Domain<BigInteger> domain, BigInteger[][] lhs, BigInteger[] rhs) {
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
                BigInteger alpha = domain.divideExact(lhs[iRow][iColumn], lhs[row][iColumn]);
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
     * Solve linear system {@code lhs.x = rhs}.
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @return solution
     * @throws ArithmeticException if the system is inconsistent or under-determined
     */
    public static BigInteger[] solve(Domain<BigInteger> domain, BigInteger[][] lhs, BigInteger[] rhs) {
        int nUnknowns = lhs[0].length;
        if (nUnknowns == 0)
            return new BigInteger[0];
        BigInteger[] result = new BigInteger[nUnknowns];
        SystemInfo info = solve(domain, lhs, rhs, result);
        if (info != SystemInfo.Consistent)
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
     * Solve linear system {@code lhs.x = rhs} and place the result to {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    public static SystemInfo solve(Domain<BigInteger> domain, BigInteger[][] lhs, BigInteger[] rhs, BigInteger[] result) {
        int nUnderDetermined = rowEchelonForm(domain, lhs, rhs);
        if (nUnderDetermined > 0)
            // under-determined system
            return SystemInfo.UnderDetermined;

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
                    return SystemInfo.Inconsistent;


        // back substitution
        for (int i = nColumns - 1; i >= 0; i--) {
            BigInteger sum = domain.getZero();
            for (int j = i + 1; j < nColumns; j++) {
                sum = domain.add(sum, domain.multiply(lhs[i][j], result[j]));
            }
            result[i] = domain.divideAndRemainder(domain.subtract(rhs[i], sum), lhs[i][i])[0];
        }
        return SystemInfo.Consistent;
    }

    /**
     * Solve linear system {@code lhs.x = rhs} and place the result to {@code result}
     * (which should be of the enough length).
     *
     * @param domain the domain
     * @param lhs    the lhs of the system
     * @param rhs    the rhs of the system
     * @param result where to place the result
     * @return system information (inconsistent, under-determined or consistent)
     */
    static SystemInfo solve(Domain<BigInteger> domain, ArrayList<BigInteger[]> lhs, ArrayList<BigInteger> rhs, BigInteger[] result) {
        return solve(domain, lhs.toArray(new BigInteger[lhs.size()][]), rhs.toArray(new BigInteger[rhs.size()]), result);
    }
}
