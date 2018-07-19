package cc.redberry.rings.scaladsl

import cc.redberry.rings.poly.multivar.GroebnerBases.HilbertSeries
import cc.redberry.rings.poly.multivar.MonomialOrder.{GREVLEX, LEX}
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

        val idealLex = ideal1.changeOrder(LEX)
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
        assert(UnivariateRing(Q, "x")("1+3*x+6*x^2+10*x^3+15*x^4+20*x^5+25*x^6+30*x^7+35*x^8+40*x^9+45*x^10+49*x^11+49*x^12+42*x^13+19*x^14") == points.hilbertSeries.numerator)

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

      val idealLex = ideal1.changeOrder(LEX)
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
      assert(UnivariateRing(Q, "x")("1+3*x+6*x^2+10*x^3+15*x^4+20*x^5+25*x^6+30*x^7+35*x^8+40*x^9+45*x^10+49*x^11+49*x^12+42*x^13+19*x^14") == points.hilbertSeries.numerator)

      val poly = x.pow(100) + y.pow(100) + z.pow(100)
      val rem = poly %% points
      assert(30 == rem.degreeSum())
    }
  }

  @Test
  def testQuotientRing1: Unit = {
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

  @Test
  def testQuotientRing2: Unit = {
    import syntax._
    val baseRing = UnivariateRing(Q, "x")
    val x = baseRing("x")

    // poly in a base ring
    val basePoly = {
      implicit val ring = baseRing
      123 * x.pow(31) + 123 * x.pow(2) + x / 2 + 1
    }

    val modulus = x.pow(2) + 1
    // poly in a quotient ring
    val quotPoly = {
      implicit val ring = AlgebraicNumberField(modulus, "t")
      123 * x.pow(31) + 123 * x.pow(2) + x / 2 + 1
    }

    assert(basePoly.degree() == 31)
    assert(quotPoly.degree() == 1)
    assert(quotPoly == basePoly % modulus)
  }

  @Test
  def testQuotientRing3: Unit = {
    import syntax._
    // base ring Q[x,y,z]
    val baseRing = MultivariateRing(Q, Array("x", "y", "z"))
    val (x, y, z) = baseRing("x", "y", "z")

    // ideal in a base ring generated by two polys <x^2 + y^12 - z, x^2*z + y^2 - 1>
    // a proper Groebner basis will be constructed automatically
    val ideal = {
      implicit val ring = baseRing
      Ideal(baseRing, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1))
    }

    // do some math in a quotient ring
    val polyQuot = {
      // quotient ring Q[x,y,z]/I
      implicit val ring = QuotientRing(baseRing, ideal)

      val poly1 = 10 * x.pow(12) + 11 * y.pow(11) + 12 * z.pow(10)
      val poly2 = x * y - y * z - z * x
      // algebraic operations performed in a quotient ring
      11 * poly1 + poly1 * poly1 * poly2
    }

    // do the same math in a base ring
    val polyBase = {
      implicit val ring = baseRing
      val poly1 = 10 * x.pow(12) + 11 * y.pow(11) + 12 * z.pow(10)
      val poly2 = x * y - y * z - z * x
      // algebraic operations performed in a base ring
      11 * poly1 + poly1 * poly1 * poly2
    }

    assert(polyQuot != polyBase)
    assert(polyQuot == polyBase %% ideal)
  }


  @Test
  def testIdeal2: Unit = {
    import syntax._

    implicit val ring = MultivariateRingZp64(17, Array("x", "y", "z"))
    val (x, y, z) = ring("x", "y", "z")

    // create ideal with two generators using GREVLEX monomial order for underlying Groebner basis
    val I = Ideal(ring, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1), GREVLEX)
    // I is proper ideal
    assert(I.isProper)

    // get computed Groebner basis
    val gb = I.groebnerBasis
    println(gb)

    // check some ideal properties
    assert(I.dimension == 1)
    assert(I.degree == 36)

    // create another ideal with only one generator
    val J = Ideal(ring, Seq(x.pow(4) * y.pow(4) + 1), GREVLEX)
    // J is principal ideal
    assert(J.isPrincipal)
    assert(J.dimension == 2)
    assert(J.degree == 8)


    val union = I union J
    // union is zero dimensional ideal
    assert(union.dimension == 0)
    // change order to LEX (elimination order)
    val eliminated = union.changeOrder(LEX)
    // system can now be solved easily
    println(eliminated)


    val intersection = I intersection J
    // intersection is still 2-dimensional
    assert(intersection.dimension == 2)
    // multiplication in this case is equal to intersection
    val times = I * J
    assert(times == intersection)


    // yet another ideal
    val K = Ideal(ring, Seq(z * x.pow(4) - z * y.pow(14) + y * z.pow(16), (x + y + z).pow(4)), GREVLEX)
    // compute complicated quotient ideal
    val quot = (I * J * K) :/ times
    assert(quot == K)
  }

  @Test
  def testIdeal3: Unit = {
    import syntax._

    // base ring in LEX order
    implicit val ring = MultivariateRing(Q, Array("x", "y", "z", "t"), LEX)
    val (x, y, z, t) = ring("x", "y", "z", "t")

    // some polynomial in a base ring order (LEX)
    val poly = x + (y ^ 2) * z + (z ^ 3) * y * t + (t ^ 4) * z * y
    assert(poly.ordering == LEX)

    // some ideal with Groebner basis computed in GREVLEX
    val idealGrevLex = Ideal(ring,
      Seq(y * x.pow(3) + z * t.pow(3) - 1,
        x * y - y * z - z * x + t.pow(3)),
      GREVLEX)
    assert(idealGrevLex.ordering == GREVLEX)

    // normal form of poly will be computed with respect to GREVLEX
    // then the result will be re-sorted according to the base ring order (LEX)
    val nfGrevLex = poly %% idealGrevLex
    assert(nfGrevLex.ordering == LEX)

    // the same ideal with Groebner basis in LEX order
    val idealLex = idealGrevLex.changeOrder(LEX)
    assert(idealLex.ordering == LEX)

    // normal form of poly will be computed with respect to LEX
    val nfLex = poly %% idealLex
    assert(nfLex.ordering == LEX)

    // Normal forms computed against LEX basis and GREVLEX basis
    // are different (although both polynomials are sorted in LEX)
    assert(nfGrevLex != nfLex)
  }

  @Test
  def testHilbertSeries1: Unit = {
    import syntax._

    // base ring in LEX order
    implicit val ring = MultivariateRing(Q, Array("x1", "x2", "x3", "x4"))
    val (x1, x2, x3, x4) = ring("x1", "x2", "x3", "x4")

    val ideal = Ideal(Seq(
      (x1 ^ 3) * (x2 ^ 3) - (x2 ^ 2) * (x4 ^ 2) - 1,
      (x1 ^ 2) * (x2 ^ 2) * (x2 ^ 3) * (x4 ^ 3) - x1 - x2 - x3 - x4,
      x1 + 2 * (x2 ^ 2) + 3 * (x3 ^ 3) + 4 * (x3 ^ 4)), GREVLEX)

    val hps: HilbertSeries = ideal.hilbertSeries

    println(hps)
    println(hps.dimension())
    //    println(hps.initialDenominatorDegree)

    println(hps.hilbertPolynomialZ())
  }

  @Test
  def testHilbertSeries2: Unit = {
    import syntax._

    implicit val ring = MultivariateRingZp64(32003, Array("x", "y", "z"))
    val (x, y, z) = ring("x", "y", "z")

    // some ideal
    val ideal = Ideal(ring, Seq(x.pow(2), y.pow(2), z.pow(2)))
    // get Hilbert-Poincare series
    val hps = ideal.hilbertSeries

    assert(hps.dimension == 0)
    assert(hps.degree == 8)

    // series numerator
    println(hps.initialNumerator)
    // reduced series numerator
    println(hps.numerator)

    // integer Hilbert polynomial
    println(hps.hilbertPolynomialZ)
    // rational Hilbert polynomial
    println(hps.hilbertPolynomial)
  }

  @Test
  def testIdeal4: Unit = {

    // ring Q[a, b, c]
    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))

    // parse some polynomials from strings
    val a = ring("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z")
    val b = ring("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4")
    val c = ring("8*x^3 + 12*y^3 + x*z^2 + 3")
    val d = ring("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3")

    // construct ideal with Groebner basis in LEX order
    val ideal = Ideal(ring, Seq(a, b, c, d), LEX)
    // it is very simple: <z^2, x, 1+4*y^3>
    println(ideal)
  }
}
