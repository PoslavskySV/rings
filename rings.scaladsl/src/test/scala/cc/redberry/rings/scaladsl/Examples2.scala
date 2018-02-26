package cc.redberry.rings.scaladsl

import cc.redberry.rings.poly.multivar.MonomialOrder.LEX
import cc.redberry.rings.primes.BigPrimes
import org.junit.Test

/**
  *
  */
class Examples2 {

  @Test
  def testIdeal1: Unit = {
    {
      def genericAssertions[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
      (implicit ring: IMultivariateRing[Term, Poly, E]): Unit = {

        import syntax._

        val (x, y, z) = ring("x", "y", "z")

        val ideal1 = Ideal(ring, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1))

        assert(1 == ideal1.dimension())
        assert(36 == ideal1.degree())
        assert(ideal1.getMonomialOrder == ring.ordering)

        val idealLex = ideal1.setMonomialOrder(LEX)
        assert(1 == idealLex.dimension())
        assert(36 == idealLex.degree())
        assert(idealLex.getMonomialOrder != ring.ordering)

        val ideal2 = Ideal(ring, Seq(x.pow(2) * y.pow(2) * z.pow(2) + x.pow(6) * y.pow(7) - 3, y.pow(2) * z.pow(2) + z.pow(6) * y.pow(7) - x))
        assert(1 == ideal2.dimension())

        val union = ideal1 ∪ ideal2
        assert(!union.isProper)
        assert(union == (ideal1 + ideal2))

        val times = ideal1 * ideal2
        assert(1 == times.dimension())

        val intersection1 = times ∩ ideal1
        assert(29 == intersection1.nBasisGenerators())

        val intersection2 = times ∩ ideal2
        assert(29 == intersection2.nBasisGenerators())

        assert(intersection1 == intersection2)

        val quot12 = ideal1 :/ ideal2
        assert(quot12 == ideal1)

        val quot21 = ideal2 :/ ideal1
        assert(7 == quot21.nBasisGenerators())
        assert(96 == quot21.degree())

        val ideal3 = Ideal(ring, Seq(x.pow(3) + y.pow(4) + z.pow(5)))
        val points = quot21 + ideal3

        assert(27 == points.nBasisGenerators())
        assert(0 == points.dimension())
        assert(UnivariateRing(Q, "x")("1+3*x+6*x^2+10*x^3+15*x^4+20*x^5+25*x^6+30*x^7+35*x^8+40*x^9+45*x^10+49*x^11+49*x^12+42*x^13+19*x^14") == points.getHilbertSeries.numerator)

        val poly = x.pow(100) + y.pow(100) + z.pow(100)
        val rem = poly %% points
        assert(30 == rem.degreeSum())
      }

      genericAssertions(MultivariateRing(Q, Array("x", "y", "z")))
      genericAssertions(MultivariateRingZp64(BigPrimes.nextPrime(1L << 59), Array("x", "y", "z")))
    }

    {
      import syntax._

      implicit val ring = MultivariateRingZp64(BigPrimes.nextPrime(1L << 59), Array("x", "y", "z"))
      val (x, y, z) = ring("x", "y", "z")

      val ideal1 = Ideal(ring, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1))

      assert(1 == ideal1.dimension())
      assert(36 == ideal1.degree())
      assert(ideal1.getMonomialOrder == ring.ordering)

      val idealLex = ideal1.setMonomialOrder(LEX)
      assert(1 == idealLex.dimension())
      assert(36 == idealLex.degree())
      assert(idealLex.getMonomialOrder != ring.ordering)

      val ideal2 = Ideal(ring, Seq(x.pow(2) * y.pow(2) * z.pow(2) + x.pow(6) * y.pow(7) - 3, y.pow(2) * z.pow(2) + z.pow(6) * y.pow(7) - x))
      assert(1 == ideal2.dimension())

      val union = ideal1 ∪ ideal2
      assert(!union.isProper)
      assert(union == (ideal1 + ideal2))

      val times = ideal1 * ideal2
      assert(1 == times.dimension())

      val intersection1 = times ∩ ideal1
      assert(29 == intersection1.nBasisGenerators())

      val intersection2 = times ∩ ideal2
      assert(29 == intersection2.nBasisGenerators())

      assert(intersection1 == intersection2)

      val quot12 = ideal1 :/ ideal2
      assert(quot12 == ideal1)

      val quot21 = ideal2 :/ ideal1
      assert(7 == quot21.nBasisGenerators())
      assert(96 == quot21.degree())


      val ideal3 = Ideal(ring, Seq(x.pow(3) + y.pow(4) + z.pow(5)))
      val points = quot21 + ideal3

      assert(27 == points.nBasisGenerators())
      assert(0 == points.dimension())
      assert(UnivariateRing(Q, "x")("1+3*x+6*x^2+10*x^3+15*x^4+20*x^5+25*x^6+30*x^7+35*x^8+40*x^9+45*x^10+49*x^11+49*x^12+42*x^13+19*x^14") == points.getHilbertSeries.numerator)

      val poly = x.pow(100) + y.pow(100) + z.pow(100)
      val rem = poly %% points
      assert(30 == rem.degreeSum())
    }
  }

  @Test
  def testQuotientRing: Unit = {
    import syntax._
    def genericAssertions[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
    (baseRing: IMultivariateRing[Term, Poly, E]): Unit = {
      val (x, y, z) = baseRing("x", "y", "z")
      val ideal = {
        implicit val __ = baseRing
        Ideal(baseRing, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1))
      }

      implicit val ring = QuotientRing(baseRing, ideal)

      val poly = x.pow(12) + y.pow(12) + z.pow(12)

      val basePoly = {
        val ring = null
        x.pow(12) + y.pow(12) + z.pow(12)
      }

      assert(poly != basePoly)
      assert(poly == basePoly %% ideal)
      assert(poly == poly %% ideal)
    }

    genericAssertions(MultivariateRing(Q, Array("x", "y", "z")))
    genericAssertions(MultivariateRingZp64(17, Array("x", "y", "z")))
  }
}
