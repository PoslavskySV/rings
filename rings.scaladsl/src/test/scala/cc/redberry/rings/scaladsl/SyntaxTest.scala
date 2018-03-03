package cc.redberry.rings.scaladsl

import cc.redberry.rings.Rational
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.PolynomialMethods.{Factor, PolynomialGCD}
import cc.redberry.rings.poly.multivar._
import cc.redberry.rings.poly.univar._
import cc.redberry.rings.primes.SmallPrimes
import org.apache.commons.math3.random.{Well1024a, Well44497a}
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.junit.Assert.{assertEquals, assertFalse, assertNotEquals, assertTrue}
import org.junit.{Assert, Ignore, Test}

import scala.language.postfixOps
import scala.reflect.ClassTag

/**
  *
  * @since 1.0
  */
class SyntaxTest {

  @Test
  def testRingOps1: Unit = {
    import syntax._
    implicit val ring = Frac(Z)
    val f = ring("2/3")
    assertEquals(ring("0"), f - f)
    assertEquals(ring("4/3"), f + f)
    assertEquals(ring("4/9"), f * f)
    assertEquals(ring("1"), f / f)
    assertEquals((ring("1"), ring("0")), f /% f)
    assertEquals(ring("0"), f % f)
    assertEquals(ring("5/3"), f ++)
    assertEquals(ring("5/3"), f + 1)
    assertEquals(ring("5/3"), 1 + f)
    assertEquals(ring("-1/3"), f --)
    assertEquals(ring("-1/3"), f - 1)
    assertEquals(ring("+1/3"), 1 - f)
    assertEquals(ring("-2/3"), -f)
    assertEquals(f, +f)
    assertEquals(f, -(-f))
    assertEquals(ring("4/3"), f * 2)
    assertEquals(ring("4/3"), 2 * f)
    assertEquals(ring("1/3"), f / 2)
    assertTrue(ring("1/3") === ring("1/3"))
    assertTrue(f =!= ring("1/3"))
    assertFalse(ring("1/3") === f)
    assertFalse(f =!= f)
    assertTrue(ring("1") === 1)
    assertFalse(ring("2") === 1)
    assertTrue(ring("1") =!= 2)
    assertFalse(ring("2") =!= 2)
  }

  @Test
  def testRingOps2: Unit = {
    def genericAssertions[E](implicit ring: Ring[E]) = {
      import syntax._
      val f = ring("2/3")
      assertEquals(ring("0"), f - f)
      assertEquals(ring("4/3"), f + f)
      assertEquals(ring("4/9"), f * f)
      assertEquals(ring("1"), f / f)
      assertEquals((ring("1"), ring("0")), f /% f)
      assertEquals(ring("0"), f % f)
      assertEquals(ring("5/3"), f ++)
      assertEquals(ring("5/3"), f + 1)
      assertEquals(ring("5/3"), 1 + f)
      assertEquals(ring("-1/3"), f --)
      assertEquals(ring("-1/3"), f - 1)
      assertEquals(ring("+1/3"), 1 - f)
      assertEquals(ring("-2/3"), -f)
      assertEquals(f, +f)
      assertEquals(f, -(-f))
      assertEquals(ring("4/3"), f * 2)
      assertEquals(ring("4/3"), 2 * f)
      assertEquals(ring("1/3"), f / 2)
      assertTrue(ring("1/3") === ring("1/3"))
      assertTrue(f =!= ring("1/3"))
      assertFalse(ring("1/3") === f)
      assertFalse(f =!= f)
      assertTrue(ring("1") === 1)
      assertFalse(ring("2") === 1)
      assertTrue(ring("1") =!= 2)
      assertFalse(ring("2") =!= 2)
    }

    genericAssertions[Rational[BigInteger]](Frac(Z))
  }

  @Test
  def testPolynomialRingOps1: Unit = {
    import syntax._
    implicit val ring = UnivariateRing(Z, "x")
    val x = ring("x")
    assertEquals(ring("0"), x - x)
    assertEquals(ring("2*x"), x + x)
    assertEquals(ring("x^2"), x * x)
    assertEquals(ring("1"), x / x)
    assertEquals((ring("1"), ring("0")), x /% x)
    assertEquals(ring("0"), x % x)
    assertEquals(ring("1 + x"), x ++)
    assertEquals(ring("1 + x"), x + 1)
    assertEquals(ring("1 + x"), 1 + x)
    assertEquals(ring("x - 1"), x --)
    assertEquals(ring("x - 1"), x - 1)
    assertEquals(ring("1 - x"), 1 - x)
    assertEquals(ring("-x"), -x)
    assertEquals(x, +x)
    assertEquals(x, -(-x))
    assertEquals(ring("2*x"), x * 2)
    assertEquals(ring("2*x"), 2 * x)
    assertEquals(ring("x"), ring("2*x") / 2)
    assertTrue(ring("x") === ring("x"))
    assertTrue(x =!= ring("x^2"))
    assertFalse(ring("2*x") === x)
    assertFalse(x =!= x)
    assertTrue(ring("(1)") === 1)
    assertFalse(ring("x") === 1)
    assertTrue(ring("x") =!= 2)
    assertFalse(ring("(2)") =!= 2)

    val y = ring("x^2")
    assertNotEquals(x, y)
    y := x
    assertEquals(x, y)

    assertEquals(ring.cfValue(0), x.cc)
    assertEquals(ring.cfValue(1), x.lc)
    assertEquals(x + 1, x + x.lc)
    assertEquals(x + 1, x.lc + x)
    assertEquals(x, x + x.cc)
    assertEquals(x, x.cc + x)
    assertEquals(x + 1, x + x.cc + 1)
    assertEquals(x + 1, x.cc + x + 1)
    assertEquals(1 + x, 1 - (x.cc - x))
    assertEquals(1 - x, 1 - (x + x.cc))
    assertEquals(1 - x, 1 - (x - x.cc))
    assertEquals(1 + x, 1 - (x.cc - x))
  }

  @Test
  def testPolynomialRingOps2: Unit = {
    def genericAssertions[Poly <: IPolynomial[Poly], E](implicit ring: IPolynomialRing[Poly, E]) = {
      import syntax._

      val x = ring("x")
      assertEquals(ring("0"), x - x)
      assertEquals(ring("2*x"), x + x)
      assertEquals(ring("x^2"), x * x)
      assertEquals(ring("1"), x / x)
      assertEquals((ring("1"), ring("0")), x /% x)
      assertEquals(ring("0"), x % x)
      assertEquals(ring("1 + x"), x ++)
      assertEquals(ring("1 + x"), x + 1)
      assertEquals(ring("1 + x"), 1 + x)
      assertEquals(ring("x - 1"), x --)
      assertEquals(ring("x - 1"), x - 1)
      assertEquals(ring("1 - x"), 1 - x)
      assertEquals(ring("-x"), -x)
      assertEquals(x, +x)
      assertEquals(x, -(-x))
      assertEquals(ring("2*x"), x * 2)
      assertEquals(ring("2*x"), 2 * x)
      assertEquals(ring("x"), ring("2*x") / 2)
      assertTrue(ring("x") === ring("x"))
      assertTrue(x =!= ring("x^2"))
      assertFalse(ring("2*x") === x)
      assertFalse(x =!= x)
      assertTrue(ring("(1)") === 1)
      assertFalse(ring("x") === 1)
      assertTrue(ring("x") =!= 2)
      assertFalse(ring("(2)") =!= 2)

      val y = ring("x^2")
      assertNotEquals(x, y)
      y := x
      assertEquals(x, y)

      assertEquals(ring.cfValue(0), x.cc)
      assertEquals(ring.cfValue(1), x.lc)
      assertEquals(x + 1, x + x.lc)
      assertEquals(x + 1, x.lc + x)
      assertEquals(x, x + x.cc)
      assertEquals(x, x * x.lc)
      assertEquals(x, x.lc * x)

      assertEquals(x, x.cc + x)
      assertEquals(x, x - x.cc)
      assertEquals(-x, x.cc - x)
      assertEquals(x.ccAsPoly, x.ccAsPoly + x.cc)
      assertEquals(x + 1, x + x.cc + 1)
      assertEquals(x + 1, x.cc + x + 1)
      assertEquals(1 + x, 1 - (x.cc - x))
      assertEquals(1 - x, 1 - (x + x.cc))
      assertEquals(1 - x, 1 - (x - x.cc))
      assertEquals(1 + x, 1 - (x.cc - x))
    }

    genericAssertions(UnivariateRing(Z, "x"))
  }

  @Test
  def testPolynomialRingOps3: Unit = {
    import syntax._

    val cfRing = UnivariateRing(Z, "z")
    implicit val ring = UnivariateRing(cfRing, "x")

    val x = ring("x")
    assertEquals(ring("0"), x - x)
    assertEquals(ring("2*x"), x + x)
    assertEquals(ring("x^2"), x * x)
    assertEquals(ring("1"), x / x)
    assertEquals((ring("1"), ring("0")), x /% x)
    assertEquals(ring("0"), x % x)
    assertEquals(ring("1 + x"), x ++)
    assertEquals(ring("1 + x"), x + 1)
    assertEquals(ring("1 + x"), 1 + x)
    assertEquals(ring("x - 1"), x --)
    assertEquals(ring("x - 1"), x - 1)
    assertEquals(ring("1 - x"), 1 - x)
    assertEquals(ring("-x"), -x)
    assertEquals(x, +x)
    assertEquals(x, -(-x))
    assertEquals(ring("2*x"), x * 2)
    assertEquals(ring("2*x"), 2 * x)
    assertEquals(ring("x"), ring("2*x") / 2)
    assertTrue(ring("x") === ring("x"))
    assertTrue(x =!= ring("x^2"))
    assertFalse(ring("2*x") === x)
    assertFalse(x =!= x)
    assertTrue(ring("(1)") === 1)
    assertFalse(ring("x") === 1)
    assertTrue(ring("x") =!= 2)
    assertFalse(ring("(2)") =!= 2)

    val y = ring("x^2")
    assertNotEquals(x, y)
    y := x
    assertEquals(x, y)

    assertEquals(ring.cfValue(0), x.cc)
    assertEquals(ring.cfValue(1), x.lc)
    assertEquals(x + 1, x + x.lc)
    assertEquals(x + 1, x.lc + x)
    assertEquals(x, x + x.cc)
    assertEquals(x, x.cc + x)
    assertEquals(x + 1, x.lc + x)
    assertEquals(x, x - x.cc)
    assertEquals(-x, x.cc + (-x))
    assertEquals(x, x * x.lc)
    assertEquals(x, (x.lc) * x)
    assertEquals(x.cc, x.cc + x.cc)

    val c1 = cfRing("1 + z")
    val c2 = cfRing("1 + z + z^2")
    val c3 = cfRing("1 + z + z^2 + z^3")
    val r = c1 + (c2 * x) + (c3 * x.pow(2))
    assertEquals(ring("(1+z)+(1+z+z^2)*x+(1+z+z^2+z^3)*x^2"), r)
  }

  @Test
  def testPolynomialRingOps4: Unit = {
    def genericAssertions[E <: IPolynomial[E]](implicit ring: IPolynomialRing[UnivariatePolynomial[E], E]) = {
      import syntax._

      val x: UnivariatePolynomial[E] = ring("x")
      assertEquals(ring("0"), x - x)
      assertEquals(ring("2*x"), x + x)
      assertEquals(ring("x^2"), x * x)
      assertEquals(ring("1"), x / x)
      assertEquals((ring("1"), ring("0")), x /% x)
      assertEquals(ring("0"), x % x)
      assertEquals(ring("1 + x"), x ++)
      assertEquals(ring("1 + x"), x + 1)
      assertEquals(ring("1 + x"), 1 + x)
      assertEquals(ring("x - 1"), x --)
      assertEquals(ring("x - 1"), x - 1)
      assertEquals(ring("1 - x"), 1 - x)
      assertEquals(ring("-x"), -x)
      assertEquals(x, +x)
      assertEquals(x, -(-x))
      assertEquals(ring("2*x"), x * 2)
      assertEquals(ring("2*x"), 2 * x)
      assertEquals(ring("x"), ring("2*x") / 2)
      assertTrue(ring("x") === ring("x"))
      assertTrue(x =!= ring("x^2"))
      assertFalse(ring("2*x") === x)
      assertFalse(x =!= x)
      assertTrue(ring("(1)") === 1)
      assertFalse(ring("x") === 1)
      assertTrue(ring("x") =!= 2)
      assertFalse(ring("(2)") =!= 2)

      val y = ring("x^2")
      assertNotEquals(x, y)
      y := x
      assertEquals(x, y)

      assertEquals(ring.cfValue(0), x.cc)
      assertEquals(ring.cfValue(1), x.lc)
      assertEquals(x + 1, x + x.lc)
      assertEquals(x + 1, x.lc + x)
      assertEquals(x, x + x.cc)
      assertEquals(x, x.cc + x)
      assertEquals(x + 1, x.lc + x)
      assertEquals(x, x - x.cc)
      assertEquals(-x, x.cc + (-x))
      assertEquals(x, x * x.lc)
      assertEquals(x, (x.lc) * x)

      assertEquals(x.cc, x.cc + x.cc)

      val cfRing: Ring[E] = ring.cfRing

      val c1 = cfRing("1 + z")
      val c2 = cfRing("1 + z + z^2")
      val c3 = cfRing("1 + z + z^2 + z^3")
      val r = c1 + (c2 * x) + (c3 * x.pow(2))
      assertEquals(ring("(1+z)+(1+z+z^2)*x+(1+z+z^2+z^3)*x^2"), r)
    }

    val cfRing = UnivariateRing(Z, "z")
    genericAssertions(UnivariateRing(cfRing, "x"))
  }

  @Test
  def testUnivariateRingOps1: Unit = {
    import syntax._

    def genericAssertions[Poly <: IUnivariatePolynomial[Poly], E <: IUnivariatePolynomial[E]](implicit ring: IUnivariateRing[Poly, E]) = {
      val x = ring("x")
      assertEquals(ring("0"), x - x)
      assertEquals(ring("2*x"), x + x)
      assertEquals(ring("x^2"), x * x)
      assertEquals(ring("1"), x / x)

      val cfRing = ring.cfRing

      assertEquals(ring cfValue 3, x.eval(3))
      assertEquals(cfRing.getZero, x.at(10))
      assertEquals(cfRing.getOne, x.at(1))

      val cfPoly = cfRing.randomElement
      // no implicit conversion
      val cfShifted = cfPoly >> 10
      assertTrue(cfShifted.degree() >= 10)
      assertTrue(cfShifted.degree() >= cfRing.perfectPowerExponent().intValue())
      assertNotEquals(cfRing valueOf cfShifted, cfShifted)

      val rPoly = ring.randomElement
      // there is implicit conversion
      val rShifted = rPoly >> 100
      assertTrue(rShifted.degree() <= ring.perfectPowerExponent().intValue())
      assertEquals(ring valueOf rShifted, rShifted)


    }

    val cfRing = GF(17, 5, "z")
    val irreducible = IrreduciblePolynomials.randomIrreduciblePolynomial(cfRing, 10, new Well44497a())
    implicit val ring = GF(irreducible, "x")
    genericAssertions(ring)
  }


  @Test
  def testGenericFunction: Unit = {
    def interpolate[Poly <: IPolynomial[Poly], E](points: Seq[E], values: Seq[E])(implicit ring: IPolynomialRing[Poly, E]) = {
      import syntax._
      points.indices
        .foldLeft(ring getZero) { case (sum, i) =>
          sum + points.indices
            .filter(_ != i)
            .foldLeft(ring getConstant values(i)) { case (product, j) =>
              implicit val cfRing = ring.cfRing
              val E: E = points(i) - points(j)
              product * (ring.`x` - points(j)) / E
            }
        }
    }

    {
      val ring = UnivariateRing(Zp(SmallPrimes.nextPrime(10000)), "x")
      val points = (1 to 50).map(BigInteger.valueOf).toArray[IntZ]
      val values: Array[IntZ] = scala.util.Random.shuffle(points.toSeq).toArray[IntZ]

      Assert.assertEquals(
        interpolate(points.toSeq, values.toSeq)(ring),
        UnivariateInterpolation.interpolateLagrange[IntZ](ring.cfRing, points, values)
      )
    }

    {

      import syntax._

      implicit val ring = UnivariateRingZp64(17, "x")
      val points = Array[Long](1, 2, 3, 5, 6, 7, 8, 9)
      val values = Array[Long](1, 2, 4, 0, 23, 1, 2, 3)

      val poly = 1 + ring("x^2 + 1").pow(7)
      assertEquals(poly, (0 to poly.degree())
        .foldLeft(poly.createZero) {
          case (sum, i) =>
            sum + poly.at(i) * ring.`x`.pow(i)
        })

      assertEquals(
        interpolate(points.toSeq, values.toSeq),
        UnivariateInterpolation.interpolateNewton(ring._cfRing, points, values)
      )
    }
  }

  final case class Ideal[
  Term <: AMonomial[Term],
  Poly <: AMultivariatePolynomial[Term, Poly],
  E](generators: Set[Poly])(implicit val ring: IMultivariateRing[Term, Poly, E]) {

    import syntax._

    private implicit val eClassTag: ClassTag[Poly] = if (generators.isEmpty) null else ClassTag(generators.head.getClass)

    private val generatorsArray: Array[Poly] = if (generators.isEmpty) null else generators.toArray[Poly]

    def isEmpty = generators.isEmpty

    def mod(poly: Poly): Poly =
      if (isEmpty)
        poly
      else
        MultivariateDivision.divideAndRemainder[Term, Poly](poly, generatorsArray: _*)(generatorsArray.length)

    def *(ideal: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] =
      if (isEmpty)
        ideal
      else
        Ideal[Term, Poly, E](for (x <- generators; y <- ideal.generators) yield x * y)

    def +(ideal: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] =
      if (isEmpty)
        ideal
      else
        Ideal[Term, Poly, E](generators ++ ideal.generators)

    def square: Ideal[Term, Poly, E] = *(this)

    def **(exponent: Int): Ideal[Term, Poly, E] = {
      if (isEmpty)
        return this

      if (exponent < 0)
        throw new IllegalArgumentException

      if (exponent == 1)
        return this

      var result: Ideal[Term, Poly, E] = Ideal.empty
      var k2p: Ideal[Term, Poly, E] = this

      var exp = exponent
      while (true) {
        if ((exp & 1) != 0)
          result = result * k2p
        exp = exp >> 1
        if (exp == 0)
          return result
        k2p = k2p * k2p
      }
      throw new RuntimeException
    }
  }

  object Ideal {
    def empty[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](implicit ring: IMultivariateRing[Term, Poly, E]) = Ideal[Term, Poly, E](Set.empty[Poly])(ring)

    def apply(generators: Seq[MultivariatePolynomialZp64])(implicit ring: MultivariateRingZp64) = Ideal[MonomialZp64, MultivariatePolynomialZp64, Long](generators.toSet)(ring)

    def apply[E](generators: Seq[MultivariatePolynomial[E]])(implicit ring: MultivariateRing[E]) = Ideal[Monomial[E], MultivariatePolynomial[E], E](generators.toSet)(ring)
  }

  class IdealOperators[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](self: Poly) {
    def mod(ideal: Ideal[Term, Poly, E]) = ideal mod self
  }

  object IdealSyntax {
    implicit def idealOperators[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](self: Poly)
    : IdealOperators[Term, Poly, E]
    = new IdealOperators[Term, Poly, E](self)

    implicit def idealOperatorsZp64(self: MultivariatePolynomialZp64)
    : IdealOperators[MonomialZp64, MultivariatePolynomialZp64, Long]
    = new IdealOperators[MonomialZp64, MultivariatePolynomialZp64, Long](self)

    implicit def idealOperatorsE[E](self: MultivariatePolynomial[E])
    : IdealOperators[Monomial[E], MultivariatePolynomial[E], E]
    = new IdealOperators[Monomial[E], MultivariatePolynomial[E], E](self)

    def idealImplicits[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
    (implicit ring: IMultivariateRing[Term, Poly, E])
    : Poly => IdealOperators[Term, Poly, E]
    = self => new IdealOperators[Term, Poly, E](self)
  }

  @Test
  def testCustomSyntax: Unit = {
    {
      import IdealSyntax._
      import syntax._

      implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))

      val ideal = Ideal(Seq(ring("x + y"), ring("y + x^2 + z"), ring("z^3")))

      val p = ring("1 + x + y + z + x*y*z + x^2 * y^2 * z^2 + 2*x*y") pow 6

      val evaled = p("x" -> ring.cfValue(2))("y" -> ring.cfValue(3))
      assertEquals(0, evaled.degree("x"))
      assertEquals(0, evaled.degree("y"))
      assertEquals(12, evaled.degree("z"))


      val mod1 = p mod ideal
      assertEquals(0, mod1.degree("x"))
      assertEquals(14, mod1.degree("y"))
      assertEquals(2, mod1.degree("z"))


      val mod2 = p mod (ideal ** 2)
      // degrees depend on order!
      //assertEquals(1, mod2.degree("x"))
      //assertEquals(16, mod2.degree("y"))
      //assertEquals(5, mod2.degree("z"))

      val mod3 = p mod (ideal + ideal ** 3)
      // degrees depend on order!
      //assertEquals(0, mod3.degree("x"))
      //assertEquals(13, mod3.degree("y"))
      //assertEquals(2, mod3.degree("z"))
    }

    def genericIdealTest[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
    (p: Poly, ring: IMultivariateRing[Term, Poly, E], ideal: Ideal[Term, Poly, E]): Unit = {
      import IdealSyntax._
      import syntax._
      implicit val pRing = ring
      implicit val _m_ = multivariateImplicits
      implicit val _i_ = idealImplicits


      val evaled = p("x" -> ring.cfValue(2))("y" -> ring.cfValue(3))
      assertEquals(0, evaled.degree("x"))
      assertEquals(0, evaled.degree("y"))
      assertEquals(12, evaled.degree("z"))


      val mod1 = p mod ideal
      assertEquals(0, mod1.degree("x"))
      assertEquals(14, mod1.degree("y"))
      assertEquals(2, mod1.degree("z"))


      val mod2 = p mod (ideal ** 2)
      // degrees depend on order!
      //assertEquals(1, mod2.degree("x"))
      //assertEquals(16, mod2.degree("y"))
      //assertEquals(5, mod2.degree("z"))

      val mod3 = p mod (ideal + ideal ** 3)
      // degrees depend on order!
      //assertEquals(0, mod3.degree("x"))
      //assertEquals(13, mod3.degree("y"))
      //assertEquals(2, mod3.degree("z"))
    }

    import syntax._
    implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))

    val ideal = Ideal(Seq(ring("x + y"), ring("y + x^2 + z"), ring("z^3")))

    val p = ring("1 + x + y + z + x*y*z + x^2 * y^2 * z^2 + 2*x*y") pow 6

    genericIdealTest(p, ring, ideal)
  }


  @Test
  def simpleUnivariateRingZp: Unit = {

    import syntax._

    implicit val ring = UnivariateRingZp64(17, "z")
    println(s"Ring $ring")

    // parse univariate polynomial from string
    val poly1 = ring("1 + 2*z^2 - z^3 + z^5 - z^17")

    // create univariate polynomial with factory method
    val poly2 = UnivariatePolynomialZ64.create(0, 1, 2, 3, 4).modulus(17)

    // another way to create poly
    val x = ring("z")
    val poly3 = 1 + (x ** 3) + 3 * x ** 9 + 9 * x ** 27

    // mix polynomials in some way
    val poly4 = (poly2 - poly1) ** 2 + (poly2 - poly1) * (poly3 - poly1) + (poly3 - poly1) ** 2

    // calculate GCD
    val gcd = PolynomialGCD(poly1 * poly3, poly2 * poly4)
    println(s"gcd : ${ring show gcd }")

    // factor some complicated poly
    val factors = Factor(poly1 * poly2 ** 2 * poly3 ** 3 * poly4 ** 4)
    println(s"factors : ${ring show factors }")
  }

  @Test
  def complicatedGaloisField: Unit = {

    val gf0 = GF(UnivariateRingZp64(17, "x")("1 + 3*x - 2*x^3 + 15*x^4 + x^5"), "x")
    val gf0Ring = UnivariateRing(gf0, "y")
    assert(IrreduciblePolynomials.irreducibleQ(gf0.getIrreducible))

    val gf1 = GF(gf0Ring("(1 + x + x^3) + (1 + x^2) * y + y^2 - 6 * y^3 + y^4 - (1 + x + x^3) * y^6 + y^7"), gf0Ring, "y")
    val gf1Ring = UnivariateRing(gf1, "z")
    assert(IrreduciblePolynomials.irreducibleQ(gf1.getIrreducible))

    val gf2 = GF(gf1Ring("((1 - x) + (1 - x^2) * y^2 + (1 - x^3) * y^3) + ((1 - x) + (1 - x^2) * y^2 + (1 - x^3) * y) * z + z^6"), gf1Ring, "z")

    implicit val ring = gf2

    println(s"ring: $ring")
    println(s"ring cardinality: ${ring.cardinality() }")
    println(s"ring characteristic: ${ring.characteristic() }")
    println(s"ring pp base: ${ring.perfectPowerBase() }")
    println(s"ring pp exponent: ${ring.perfectPowerExponent() }")

    val rndPoly = ring.randomElement(new Well1024a())
    println(ring show rndPoly)
    println(ring show Factor(rndPoly))
  }

  @Test
  def monomials: Unit = {
    implicit val ring = MultivariateRingZp64(3, Array("x", "y", "z"))

    val m1 = MonomialZp64(3, "x" -> 1, "y" -> 2)
    val m2 = MonomialZp64(3, "x" -> 1, "y" -> 3)
  }

  @Test
  def integers1: Unit = {
    import syntax._

    val ring = Zp(2)
    val a = ring("2")
    val b = ring("3")
    assertEquals(ring(5), a + b)

    implicit val implicitRing = ring
    assertEquals(ring(1), a + b)
  }

  @Test
  def testUnivariateRingOps2: Unit = {
    import syntax._

    val ring = UnivariateRingZp64(17, "x")

    val divider = ring("1 + x + x^2 + x^3").pow(15)
    val list = (1 until 1000)
      .map(_ => ring.randomElement(divider.degree() * 3 / 2, divider.degree() * 2, new Well1024a()))

    // precomputed Newton inverses
    implicit val invMod: PrecomputedInverse[UnivariatePolynomialZp64] = divider.precomputedInverses

    var (timePlain, timeFast) = (new DescriptiveStatistics(), new DescriptiveStatistics())

    var counter = 0
    for (el <- list) {
      counter += 1
      if (counter == list.length / 10) {
        timePlain.clear()
        timeFast.clear()
      }

      var start = System.nanoTime()
      // quotient and remainder using built-in methods
      val divRemPlain = el /% divider
      timePlain.addValue(System.nanoTime() - start)


      start = System.nanoTime()
      // quotient and remainder computed using fast
      // algorithm with Newton iterations
      val divRemFast = el /%% divider
      timeFast.addValue(System.nanoTime() - start)

      assert(divRemPlain == divRemFast)
    }

    println(timePlain)
    println(timeFast)
  }

  @Test
  def testMultivariateRingOps2: Unit = {
    import syntax._
    implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))

    val poly = ring("x + y + x")
    val p1 = ring("x - y + x")
    val p2 = ring("x + y - x")

    println(poly /%/% (p1, p2))
    println(poly /%/% (p1, p2, p2))
    println(poly /%/% (p1, p2, p2, p1))
  }

  @Ignore //issue #23
  @Test
  def testParse1: Unit = {

    implicit val gf = GF(17, 3, "t")
    implicit val ring = MultivariateRing(gf, Array("x", "y", "z"))
    val t = ring("t")
  }
}
