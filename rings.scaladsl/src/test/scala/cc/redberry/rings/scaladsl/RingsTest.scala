package cc.redberry.rings.scaladsl

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.IPolynomial
import cc.redberry.rings.poly.PolynomialMethods.{Factor, PolynomialGCD}
import cc.redberry.rings.poly.multivar._
import cc.redberry.rings.poly.univar.{IrreduciblePolynomials, UnivariateGCD, UnivariateInterpolation, UnivariatePolynomialZ64}
import cc.redberry.rings.primes.SmallPrimes
import cc.redberry.rings.scaladsl.Rings._
import org.apache.commons.math3.random.Well1024a
import org.junit.{Assert, Test}

import scala.reflect.ClassTag

/**
  *
  *
  * @since 1.0
  */
class RingsTest {

  import implicits._

  @Test
  def simpleUnivariateRingZp: Unit = {

    implicit val ring = UnivariateRingZp64(17, "z")
    println(s"Ring $ring")

    // parse univariate polynomial from string
    val poly1 = ring("1 + 2*z^2 - z^3 + z^5 - z^17")

    // create univariate polynomial with factory method
    val poly2 = UnivariatePolynomialZ64.create(0, 1, 2, 3, 4).modulus(17)

    // another way to create poly
    val x = parse("z")
    val poly3 = 1 + (x ** 3) + 3 * x ** 9 + 9 * x ** 27

    // mix polynomials in some way
    val poly4 = (poly2 - poly1) ** 2 + (poly2 - poly1) * (poly3 - poly1) + (poly3 - poly1) ** 2

    // calculate GCD
    val gcd = PolynomialGCD(poly1 * poly3, poly2 * poly4)
    println(s"gcd : ${ring show gcd}")

    // factor some complicated poly
    val factors = Factor(poly1 * poly2 ** 2 * poly3 ** 3 * poly4 ** 4)
    println(s"factors : ${ring show factors}")
  }


  @Test
  def simpleMultivariateRingZ: Unit = {
    import cc.redberry.rings.poly.PolynomialMethods.{Factor, PolynomialGCD}
    import implicits._

    implicit val ring = MultivariateRing(Rings.Z, Array("x", "y", "z"))
    println(s"Ring $ring")

    val poly1 = ring("x + y + z") ** 10 + 1
    val poly2 = ring("x - y^2 + z^3") ** 10 + 1

    // calculate GCD
    val gcd = PolynomialGCD(poly1, poly2)
    println(s"gcd: ${ring show gcd}")

    // factor some complicated poly
    val factors = Factor(poly1 * poly2)
    println(s"Factor decomposition: ${ring show factors}")
  }

  @Test
  def complicatedGaloisField: Unit = {
    import implicits._

    val gf0 = GF(UnivariateRingZp64(17, "x")("1 + 3*x - 2*x^3 + 15*x^4 + x^5"), "x")
    val gf0Ring = UnivariateRing(gf0, "y")
    assert(IrreduciblePolynomials.irreducibleQ(gf0.getIrreducible))

    val gf1 = GF(gf0Ring("(1 + x + x^3) + (1 + x^2) * y + y^2 - 6 * y^3 + y^4 - (1 + x + x^3) * y^6 + y^7"), gf0Ring, "y")
    val gf1Ring = UnivariateRing(gf1, "z")
    assert(IrreduciblePolynomials.irreducibleQ(gf1.getIrreducible))

    val gf2 = GF(gf1Ring("((1 - x) + (1 - x^2) * y^2 + (1 - x^3) * y^3) + ((1 - x) + (1 - x^2) * y^2 + (1 - x^3) * y) * z + z^6"), gf1Ring, "z")

    implicit
    val ring = gf2

    println(s"ring: $ring")
    println(s"ring cardinality: ${ring.cardinality()}")
    println(s"ring characteristic: ${ring.characteristic()}")
    println(s"ring pp base: ${ring.perfectPowerBase()}")
    println(s"ring pp exponent: ${ring.perfectPowerExponent()}")

    val rndPoly = ring.randomElement(new Well1024a())
    println(rndPoly)
    println(Factor(rndPoly))
  }

  @Test
  def complicatedMultivariateRing: Unit = {
    import implicits._

    val coefficientRing = GF(UnivariateRingZp64(41, "x")("1 + x + x^2"), "z")

    implicit val ring = MultivariateRing(coefficientRing, Array("a", "b", "c"), MonomialOrder.LEX)

    println(s"Coefficient ring: $coefficientRing")

    // parse univariate polynomial from string
    val poly1 = parse("(1 + x + x^2)*a*b*c + a^2 + (x)*c^3 + (13 + 12*x + 4*x^6)*b^2*c^3 + (2 + x)")
    val poly2 = parse("(1 - x - x^2)*a*b*c - a^2 - (x)*c^3 + (13 - 12*x + 4*x^6)*b^2*c^3 + (2 + x^87)")
    val poly = poly1 * poly2 ** 2 * (poly1 + poly2) ** 3

    val f = coefficientRing("x") + coefficientRing("2*x")
    println(f)

    println(Factor(poly))
  }


  @Test
  def interpolation1: Unit = {
    implicit val ring = UnivariateRing(Rings.Zp(SmallPrimes.nextPrime(10000)), "x")

    val points = (1 to 50).map(BigInteger.valueOf).toArray[Integer]
    val values: Array[Integer] = scala.util.Random.shuffle(points.toSeq).toArray[Integer]

    Assert.assertEquals(
      interpolate(points.toSeq, values.toSeq),
      UnivariateInterpolation.interpolateLagrange[Integer](ring.coefficientDomain, points, values)
    )
  }

  @Test
  def interpolation2: Unit = {
    import implicits._
    implicit val ring = UnivariateRingZp64(17, "x")
    val points = Array[Long](1, 2, 3, 5, 6, 7, 8, 9)
    val values = Array[Long](1, 2, 4, 0, 23, 1, 2, 3)

    val poly = 1 + ring("x^2 + 1") ** 7
    println(poly == (0 to poly.degree()).foldLeft(poly createZero()) { case (sum, i) => sum + poly(i) * (ring.`x` ** i) })

    println(poly.cc())
    println(interpolate(points.toSeq, values.toSeq))
    println(UnivariateInterpolation.interpolateNewton(ring.coefficientDomain, points, values))
  }

  def interpolate[Poly <: IPolynomial[Poly], E](points: Seq[E], values: Seq[E])(implicit ring: PolynomialRing[Poly, E]) = {
    import implicits.generic.polynomials._
    points.indices
      .foldLeft(ring getZero) { case (sum, i) =>
        sum + points.indices
          .filter(_ != i)
          .foldLeft(ring getConstant values(i)) { case (product, j) =>
            product * (ring.`x` - points(j)) / (points(i) - points(j))
          }
      }
  }


  @Test
  def ideal1: Unit = {
    import implicits._

    implicit val ring = MultivariateRing(Rings.Z, Array("x", "y", "z"))

    val ideal = Ideal(Seq(ring("x + y"), ring("y + x^2 + z"), ring("z^3")))

    val p = ring("1 + x + y + z + x*y*z + x^2 * y^2 * z^2 + 2*x*y") ** 6

    println(p("x" -> 2)("y" -> 3))
    //    println(p mod ideal)
    //    println(p mod ideal ** 2)
    //    println(p mod (ideal + ideal ** 3))
    //    println(p mod ideal ** 4)
  }

  case class IdealOperators[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](self: Poly) {
    def mod(ideal: Ideal[Term, Poly, E]) = ideal mod self
  }

  implicit def withIdealOperators[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](self: Poly)
  : IdealOperators[Term, Poly, E]
  = IdealOperators[Term, Poly, E](self)

  implicit def withIdealOperatorsZp64(self: MultivariatePolynomialZp64)
  : IdealOperators[MonomialZp64, MultivariatePolynomialZp64, Long]
  = IdealOperators[MonomialZp64, MultivariatePolynomialZp64, Long](self)

  case class Ideal[
  Term <: DegreeVector[Term],
  Poly <: AMultivariatePolynomial[Term, Poly],
  E] private(generators: Set[Poly])(implicit val ring: IMultivariateRing[Term, Poly, E]) {

    import implicits.generic.polynomials._

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
    def empty[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](implicit ring: IMultivariateRing[Term, Poly, E]) = Ideal[Term, Poly, E](Set.empty[Poly])(ring)

    def apply(generators: Seq[MultivariatePolynomialZp64])(implicit ring: MultivariateRingZp64) = Ideal[MonomialZp64, MultivariatePolynomialZp64, Long](generators.toSet)(ring)

    def apply[E](generators: Seq[MultivariatePolynomial[E]])(implicit ring: MultivariateRing[E]) = Ideal[Monomial[E], MultivariatePolynomial[E], E](generators.toSet)(ring)
  }

  @Test
  def xxxx: Unit = {
    import implicits._

    val polyRing = MultivariateRingZp64(2, Array("x", "y", "z"))

    implicit val ring = Rationals(polyRing)

    val a = ring("(x + y + z)/(x - y)")
    val b = ring("(x^2 - y + z)/(x + y)")
    val c = ring("(x^2 + y - z^2)/(x - y)")
    val sum = a * b ** 2 / c - a * c / b ** 3 + b ** 4 * c / a

    //        println(sum)

    assert(sum % a === 0)
    println(Factor(sum.denominator))
    //    println(Factor(sum.numerator))
  }

  @Test
  def zxzxjchzkcj: Unit = {
    import implicits._

    implicit val ring = MultivariateRingZp64(30203, Array("a", "b", "c"))

    println(Factor(ring("20407*c+15337*b^2*c+15309*a*b^3+21682*a*b^3*c^2+16926*a^3*b^3*c")))
    println(Factor(ring("1398*b^2*c^2+28996*a*c^3+8241*a*b^2*c+14631*a^2*b*c+19151*a^3*b^3")))
    println(Factor(ring("3422*b^2*c^2+23893*a+20695*a*b*c+18374*a*b^2+29283*a^3*b^2*c")))
    println(Factor(ring("7620*c^2+24026*b*c^2+6319*a^2")))
  }
}
