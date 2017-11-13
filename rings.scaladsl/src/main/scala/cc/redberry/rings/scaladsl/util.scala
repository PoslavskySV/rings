package cc.redberry.rings.scaladsl

/**
  * @since 1.0
  */
object util {

  /**
    * Evaluates fun and measures the time elapsed (in seconds)
    *
    * @return (time, result)
    */
  def timing[R](fun: => R): (Double, R) = {
    val start = System.nanoTime()
    val r: R = fun
    val timing = (System.nanoTime() - start).toDouble / 1000000000.0
    (timing, r)
  }
}
