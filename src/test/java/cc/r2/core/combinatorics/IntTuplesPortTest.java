package cc.r2.core.combinatorics;

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IntTuplesPortTest {

    @Test
    public void test1() {
        IntTuplesPort port = new IntTuplesPort(3, 3, 3);
        int count = 0;
        while (port.take() != null)
            ++count;
        Assert.assertEquals(count, 27);
    }

    @Test
    public void test2() {
        IntTuplesPort port = new IntTuplesPort(4, 4, 4, 4);
        int count = 0;
        while (port.take() != null)
            ++count;
        Assert.assertEquals(count, 256);
    }

    @Test
    public void testLUD() throws Exception {
        int[][] results = {
                {0, 0, 0},
                {0, 0, 1},
                {0, 1, 0},
                {0, 1, 1},
                {0, 2, 0},
                {0, 2, 1},
                {1, 0, 0},
                {1, 0, 1},
                {1, 1, 0},
                {1, 1, 1},
                {1, 2, 0},
                {1, 2, 1}
        };

        int[] luds = {0, 2, 1, 2, 1, 2, 0, 2, 1, 2, 1, 2};

        IntTuplesPort port = new IntTuplesPort(2, 3, 2);
        int i = 0;
        int[] r;
        while ((r = port.take()) != null) {
            Assert.assertArrayEquals(results[i], r);
            Assert.assertEquals(luds[i++], port.getLastUpdateDepth());
        }
    }

    @Test
    public void test3() {
        IntTuplesPort port = new IntTuplesPort(3, 2, 2);
        int[] c;
        while ((c = port.take()) != null) {
            System.out.println(Arrays.toString(c));
        }
    }
}
