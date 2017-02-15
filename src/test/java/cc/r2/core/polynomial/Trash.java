package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

/**
 * Created by poslavsky on 14/02/2017.
 */
public class Trash {
    public static abstract class AbstractClass<T extends AbstractClass> {
        long[] data;

        public AbstractClass(long[] data) {
            this.data = data;
        }

        private final T self = (T) this;

        T ensureCapacity(int newLen) {
            if (data.length < newLen) {
                int oldLen = data.length;
                data = Arrays.copyOf(data, newLen);
                Arrays.fill(data, oldLen, newLen, 0);
            }
            return self;
        }

        public abstract T add(T arr1);

        public abstract T create(long[] data);
    }

    public static final class AImpl1 extends AbstractClass<AImpl1> {
        public AImpl1(long[] data) {
            super(data);
        }

        @Override
        public AImpl1 create(long[] data) {
            return new AImpl1(data);
        }

        @Override
        public AImpl1 add(AImpl1 arr1) {
            ensureCapacity(arr1.data.length);
            for (int i = 0; i < arr1.data.length; i++)
                data[i] += arr1.data[i];
            return this;
        }
    }

    public static final class AImpl2 extends AbstractClass<AImpl2> {
        public AImpl2(long[] data) {
            super(data);
        }

        @Override
        public AImpl2 create(long[] data) {
            return new AImpl2(data);
        }

        @Override
        public AImpl2 add(AImpl2 arr1) {
            ensureCapacity(arr1.data.length);
            for (int i = 0; i < arr1.data.length; i++)
                data[i] -= arr1.data[i] / 2;
            return this;
        }
    }

    public static final class BImpl1 {
        long[] data;

        public BImpl1(long[] data) {
            this.data = data;
        }

        public BImpl1 add(BImpl1 arr1) {
            if (data.length < arr1.data.length)
                data = Arrays.copyOf(data, arr1.data.length);
            for (int i = 0; i < arr1.data.length; i++)
                data[i] += arr1.data[i];
            return this;
        }
    }

    public static final class BImpl2 {
        long[] data;

        public BImpl2(long[] data) {
            this.data = data;
        }

        public BImpl2 add(BImpl2 arr1) {
            if (data.length < arr1.data.length)
                data = Arrays.copyOf(data, arr1.data.length);
            for (int i = 0; i < arr1.data.length; i++)
                data[i] -= arr1.data[i] / 2;
            return this;
        }
    }

    @SuppressWarnings("unchecked")
    public static <T extends AbstractClass<T>> T aJob(long[] data, T doer) {
        T r = doer.create(data);
        r = r.add(doer);
        for (int i = 0; i < data.length; i++) {
            r = r.add(r).add(doer);
            ++r.data[i];
        }
        return r;
    }

    public static BImpl1 bJob1(long[] data, BImpl1 doer) {
        BImpl1 r = new BImpl1(data);
        r = r.add(doer);
        for (int i = 0; i < data.length; i++) {
            r = r.add(r).add(doer);
            ++r.data[i];
        }
        return r;
    }

    public static BImpl2 bJob2(long[] data, BImpl2 doer) {
        BImpl2 r = new BImpl2(data);
        r = r.add(doer);
        for (int i = 0; i < data.length; i++) {
            r = r.add(r).add(doer);
            ++r.data[i];
        }
        return r;
    }


    @Test
    public void benchmarkReflection() throws Exception {
        DescriptiveStatistics aImpl1 = new DescriptiveStatistics(),
                aImpl2 = new DescriptiveStatistics(),
                bImpl1 = new DescriptiveStatistics(),
                bImpl2 = new DescriptiveStatistics();

        List<DescriptiveStatistics> stats = Arrays.asList(aImpl1, aImpl2, bImpl1, bImpl1);
        RandomGenerator rnd = new Well1024a();
        int len = 50;
        for (int i = 0; i < 120_000; i++) {
            if (i == 10_000)
                stats.forEach(DescriptiveStatistics::clear);

            long[] data = new long[len];
            long[] baseData = new long[len];
            for (int j = 0; j < len; j++) {
                data[j] = rnd.nextLong();
                baseData[j] = rnd.nextLong();
            }

            BImpl1 b1 = new BImpl1(baseData.clone());
            BImpl2 b2 = new BImpl2(baseData.clone());
            AImpl1 a1 = new AImpl1(baseData.clone());
            AImpl2 a2 = new AImpl2(baseData.clone());

            long start;

            start = System.nanoTime();
            long[] r_b1 = bJob1(data, b1).data;
            bImpl1.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            long[] r_b2 = bJob2(data, b2).data;
            bImpl2.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            long[] r_a1 = aJob(data, a1).data;
            aImpl1.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            long[] r_a2 = aJob(data, a2).data;
            aImpl2.addValue(System.nanoTime() - start);

            assertArrayEquals(r_a1, r_b1);
            assertArrayEquals(r_a2, r_b2);
        }

        System.out.println("aImpl1: " + aImpl1.getPercentile(50));
        System.out.println("aImpl2: " + aImpl2.getPercentile(50));
        System.out.println("bImpl1: " + bImpl1.getPercentile(50));
        System.out.println("bImpl2: " + bImpl2.getPercentile(50));
    }
}
