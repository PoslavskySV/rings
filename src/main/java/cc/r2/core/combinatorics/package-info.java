/**
 * <p>Provides a number of useful combinatorial algorithms. It includes two main types of classes:
 * classes aimed on the enumeration of a particular kinds of combinations, and classes providing facilities to
 * operate with single permutations and symmetric groups.</p>
 *
 * <p><b><i>Enumerating combinations.</i></b> All of the classes, which implements the algorithms of combinations enumeration
 * are strictly follows a common pattern of output port, which is specified in
 * {@link cc.r2.core.combinatorics.IntCombinatorialPort}. The calculation of the next combination occurs
 * strictly on the invocation of method {@link cc.r2.core.combinatorics.IntCombinatorialPort#take()} and the
 * returned array is always the same reference. Some of these classes are also implements {@link java.util.Iterator} and {@link java.lang.Iterable} interfaces for convenience.
 * <table>
 *     <tr>
 *         <td><b>List of enumeration algorithms:</b></td>
 *         <td></td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntPermutationsGenerator}</td>
 *         <td>Enumerates all permutations of dimension N ( N! permutations). </td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntCombinationsGenerator}</td>
 *         <td>Enumerates all combinations of K elements chosen as N ( N!/(K!(N-K)!) combinations). </td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntCombinationPermutationGenerator}</td>
 *         <td>Enumerates all combinations with permutations of K elements chosen as N ( N!/(N-K)! combinations). </td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntDistinctTuplesPort}</td>
 *         <td>Enumerates all distinct N-tuples, which can be chosen from {@code N} sets of integers. </td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntTuplesPort}</td>
 *         <td>Enumerates all N-tuples, which can be chosen from {@code N} sets of integers of the form
 *               <i>array</i><sub>i</sub> = [0, 1, 2, ..., K<sub>i</sub>]. </td>
 *     </tr>
 *     <tr>
 *         <td>{@link cc.r2.core.combinatorics.IntPriorityPermutationsGenerator}</td>
 *         <td>Enumerates all permutations of dimension N ( N! permutations) and allows to affect on the
 *              enumeration order.</td>
 *     </tr>
 * </table>
 * </p>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
package cc.r2.core.combinatorics;
