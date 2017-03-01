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
package cc.r2.core.util;

import java.util.Iterator;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public interface OutputPort<T> {
    T take();

    final class Singleton<T> implements OutputPort<T> {
        private T element;

        public Singleton(T element) {
            this.element = element;
        }

        @Override
        public T take() {
            T newElement = element;
            element = null;
            return newElement;
        }
    }

    final class PortIterator<T> implements Iterator<T> {
        private final OutputPort<T> opu;
        private T next;

        public PortIterator(OutputPort<T> opu) {
            this.opu = opu;
        }

        @Override
        public boolean hasNext() {
            return (next = opu.take()) != null;
        }

        @Override
        public T next() {
            return next;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Not supported yet.");
        }
    }

    final class PortIterable<T> implements Iterable<T> {
        private final OutputPort<T> opu;

        public PortIterable(OutputPort<T> opu) {
            this.opu = opu;
        }

        @Override
        public Iterator<T> iterator() {
            return new PortIterator<>(opu);
        }
    }
}
