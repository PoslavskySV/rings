package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
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

    static String[] defaultVars(int nVars) {
        char v = 'a';
        String[] vars = new String[nVars];
        for (int i = 0; i < nVars; i++)
            vars[i] = Character.toString(v++);
        return vars;
    }
}
