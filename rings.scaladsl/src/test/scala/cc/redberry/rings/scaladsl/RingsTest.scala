package cc.redberry.rings.scaladsl

import cc.redberry.rings.scaladsl.syntax._
import org.junit.Assert.assertEquals
import org.junit.Test

/**
  *
  */
class RingsTest {
  @Test
  def simpleRings: Unit = {
    {
      implicit val ring = Z
      assertEquals(ring(3), ring("1 + 4 / 2"))
      assertEquals(ring(3), ring(1) + ring(4) / ring(2))
      val (x, y) = ring(3, 17)
      assertEquals(x - y, -ring("14"))
    }

    {
      implicit val ring = Q
      assertEquals(ring(3) / 2, ring("(1 + 4 / 2) / 2"))
      assertEquals(ring(3) / ring(4), (ring(1) + ring(4) / ring(2)) / 4)
      val (x, y) = ring(3, 17)
      assertEquals(x - y, -ring("14"))
    }

    {
      implicit val ring = Zp(17)
      assertEquals(ring(3) / 2, ring("(1 + 4 / 2) / 2"))
      assertEquals(ring(3) / ring(4), (ring(1) + ring(4) / ring(2)) / 4)
      val (x, y) = ring(3, 17)
      assertEquals(x - y, -ring("14"))
    }
  }

  @Test
  def GaloisFields64: Unit = {
    implicit val gf = GF(UnivariateRingZp64(17, "t")("13 + 8*t + t^3"), "t")
    assertEquals("(Z/17)[t]/<13+8*t+t^3>", gf.toString)
    val t = gf("t")
    assertEquals(1 + t + t.pow(3) / (t - 1 - t.pow(66)), gf("1 + t + t^3 / (t - 1 - t^66)"))
  }

  @Test
  def GaloisFieldsE: Unit = {
    implicit val cfRing = GF(17, 3, "t")

    val uRing = UnivariateRing(cfRing, "x")
    val irreducible = uRing("4*t+11*t^2+(7*t+15*t^2)*x+(5+7*t+2*t^2)*x^2+x^3")

    implicit val gf = GF(irreducible, "x")
    assertEquals("((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x]/<4*t+11*t^2+(7*t+15*t^2)*x+(5+7*t+2*t^2)*x^2+x^3>", gf.toString)

    val el = gf("t + x")
    assertEquals("t+x", gf stringify el)

    val (x, t) = gf("x", "t")
    assertEquals(t + x - x.pow(2) / t - x / (t + x), gf("t + x  - x^2 / t - x / (t + x)"))
  }

  @Test
  def MultivarRings: Unit = {
    implicit val gf = GF(UnivariateRingZp64(17, "t")("13 + 8*t + t^3"), "t")
    implicit val polyRing = MultivariateRing(gf, Array("x", "y", "z"))
    implicit val fracRing = Frac(polyRing)
    implicit val wRing = UnivariateRing(fracRing, "W")

    val poly = wRing("(t / x + W - W^2 / (x - y) + z) * (W - 2)")
    val factors = wRing factor poly
    assert(wRing(wRing stringify factors) == poly)
    println(wRing)
  }
}
