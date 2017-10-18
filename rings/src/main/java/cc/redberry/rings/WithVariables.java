package cc.redberry.rings;

/**
 * @since 1.0
 */
public interface WithVariables {
    /**
     * Returns string representation of this using specified string representation for variables.
     *
     * @param variables string values of variables
     * @return string representation of this
     */
    String toString(String[] variables);

    /**
     * Sequence of strings "a", "b", "c" etc.
     *
     * @param nVars number of variable
     */
    static String[] defaultVars(int nVars) {
        if (nVars == 1)
            return new String[]{"x"};
        if (nVars == 2)
            return new String[]{"x", "y"};
        if (nVars == 3)
            return new String[]{"x", "y", "z"};

        String[] vars = new String[nVars];
        for (int i = 1; i <= nVars; i++)
            vars[i - 1] = "x" + i;
        return vars;
    }
}
