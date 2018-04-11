package cc.redberry.rings.scaladsl

import org.junit.Test

/**
  *
  * @since 2.3
  */
class ConversionsTest {
  @Test
  def test1: Unit = {
    import Conversions._

    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))
    val poly = ring("x + y + z")

    println(split(ring, 1, 2))
    println(fromUnivariate(asUnivariate(ring, 1), 1))
  }
}
