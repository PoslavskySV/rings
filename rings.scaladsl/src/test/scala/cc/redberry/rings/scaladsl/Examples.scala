package cc.redberry.rings.scaladsl

import org.junit.Test

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
class Examples {
  @Test
  def test1: Unit = {
    implicit val ring = GF(3, 4, "x")

    println(UnivariatePolynomial(1, 2, 3))

    println(UnivariatePolynomial(ring("x"), ring("1 + x")))
  }
}
