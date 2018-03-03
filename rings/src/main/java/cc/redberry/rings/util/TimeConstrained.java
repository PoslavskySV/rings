package cc.redberry.rings.util;

import java.util.concurrent.*;

/**
 *
 */
public class TimeConstrained {
    /**
     * Runs lambda, stopping after specified number of milliseconds
     *
     * @param lambda   code block
     * @param millis   time constrained in milliseconds
     * @param failExpr returns this if the time constraint is not met
     */
    public static <T> T timeConstrained0(Callable<T> lambda, long millis, T failExpr)
            throws InterruptedException, ExecutionException {
        ExecutorService executor = Executors.newSingleThreadExecutor();
        FutureTask<T> task = new FutureTask<>(lambda);
        executor.execute(task);
        try {
            return task.get(millis, TimeUnit.MILLISECONDS);
        } catch (TimeoutException e) {
            return failExpr;
        } finally {
            task.cancel(true);
            executor.shutdown();
        }
    }

    /**
     * Runs lambda, stopping after specified number of milliseconds
     *
     * @param lambda   code block
     * @param millis   time constrained in milliseconds
     * @param failExpr returns this if the time constraint is not met
     */
    public static <T> T timeConstrained(Callable<T> lambda, long millis, T failExpr) {
        try {
            return timeConstrained0(lambda, millis, failExpr);
        } catch (InterruptedException | ExecutionException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Runs lambda, stopping after specified number of milliseconds
     *
     * @param lambda code block
     * @param millis time constrained in milliseconds
     */
    public static <T> T timeConstrained(Callable<T> lambda, long millis)
            throws InterruptedException, ExecutionException {
        return timeConstrained0(lambda, millis, null);
    }
}
