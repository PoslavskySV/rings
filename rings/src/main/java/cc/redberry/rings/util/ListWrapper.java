package cc.redberry.rings.util;

import java.util.*;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;
import java.util.stream.Stream;

/** A simple list wrapper */
public class ListWrapper<Poly>
        extends AbstractList<Poly> {
    /** Inner list */
    public final List<Poly> list;

    public ListWrapper(List<Poly> list) { this.list = list; }

    @Override
    public boolean isEmpty() {
        return list.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return list.contains(o);
    }

    @Override
    public Object[] toArray() {
        return list.toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return list.toArray(a);
    }

    @Override
    public boolean remove(Object o) {
        return list.remove(o);
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return list.containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends Poly> c) {
        return list.addAll(c);
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        return list.removeAll(c);
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        return list.retainAll(c);
    }

    @Override
    public void replaceAll(UnaryOperator<Poly> operator) {
        list.replaceAll(operator);
    }

    @Override
    public void sort(Comparator<? super Poly> c) {
        list.sort(c);
    }

    @Override
    public Spliterator<Poly> spliterator() {
        return list.spliterator();
    }

    @Override
    public boolean removeIf(Predicate<? super Poly> filter) {
        return list.removeIf(filter);
    }

    @Override
    public Stream<Poly> stream() {
        return list.stream();
    }

    @Override
    public Stream<Poly> parallelStream() {
        return list.parallelStream();
    }

    @Override
    public void forEach(Consumer<? super Poly> action) {
        list.forEach(action);
    }

    @Override
    public boolean add(Poly poly) {
        return list.add(poly);
    }

    @Override
    public Poly set(int index, Poly element) {
        return list.set(index, element);
    }

    @Override
    public void add(int index, Poly element) {
        list.add(index, element);
    }

    @Override
    public Poly remove(int index) {
        return list.remove(index);
    }

    @Override
    public int indexOf(Object o) {
        return list.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return list.lastIndexOf(o);
    }

    @Override
    public void clear() {
        list.clear();
    }

    @Override
    public boolean addAll(int index, Collection<? extends Poly> c) {
        return list.addAll(index, c);
    }

    @Override
    public Iterator<Poly> iterator() {
        return list.iterator();
    }

    @Override
    public ListIterator<Poly> listIterator() {
        return list.listIterator();
    }

    @Override
    public ListIterator<Poly> listIterator(int index) {
        return list.listIterator(index);
    }

    @Override
    public List<Poly> subList(int fromIndex, int toIndex) {
        return list.subList(fromIndex, toIndex);
    }

    @Override
    public boolean equals(Object o) {
        return list.equals(o);
    }

    @Override
    public int hashCode() {
        return list.hashCode();
    }

    @Override
    protected void removeRange(int fromIndex, int toIndex) {
        super.removeRange(fromIndex, toIndex);
    }

    @Override
    public Poly get(int index) {
        return list.get(index);
    }

    @Override
    public int size() {
        return list.size();
    }

    @Override
    public String toString() {
        return list.toString();
    }
}