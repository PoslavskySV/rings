/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.rings.util;

import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class IntTimSortTest {
    private static final int K = 200;

    @Test
    public void testSomeMethod() {
        int pivot = K / 6;
        Random r = new Random();
        for (int n = 0; n < 1000; ++n) {
            int[] main = new int[K];
            int[] perm = new int[K];
            int counter = 0;
            for (int i = 0; i < K; ++i) {
                perm[i] = main[i] = r.nextInt(K / 4);
                if (main[i] == pivot)
                    perm[i] = counter++;
            }
            IntTimSort.sort(main, perm);
            counter = 0;
            for (int i = 0; i < K; ++i) {
                if (i < K - 1)
                    assertTrue(main[i] <= main[i + 1]);
                if (main[i] != pivot)
                    assertEquals(main[i], perm[i]);
                else
                    assertEquals(counter++, perm[i]);
            }
        }
    }
}
