package cc.redberry.rings.scaladsl

import cc.redberry.rings
import cc.redberry.rings.poly.multivar.MultivariatePolynomial
import org.junit.Test

import scala.language.postfixOps
import scala.util.Random

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
class Examples {
  @Test
  def test1: Unit = {

    import syntax._
    // when parsing "x" will be considered as the "first variable"
    // and "y" as "the second" => in the result the particular
    // names "x" and "y" are erased
    val poly1 = MultivariatePolynomial.parse("x^2 + x*y", "x", "y")
    // parse the same polynomial but using "a" and "b" instead of "x" and "y"
    val poly2 = MultivariatePolynomial.parse("a^2 + a*b", "a", "b")
    // polynomials are equal (no matter which variable names were used when parsing)
    assert(poly1 == poly2)
    // degree in the first variable
    assert(poly1.degree(0) == 2)
    // degree in the second variable
    assert(poly1.degree(1) == 1)

    // this poly differs from poly2 since now "a" is "the second"
    // variable and "b" is "the first"
    val poly3 = MultivariatePolynomial.parse("a^2 + a*b", "b", "a")
    assert(poly3 != poly2)
    // swap the first and the second variables and the result is equal to poly2
    assert(poly3.swapVariables(0, 1) == poly2)


    // the default toString() will use the default
    // variables "a", "b", "c"  and so on (alphabetical)
    // the result will be "a*b + a^2"
    println(poly1)
    // specify which variable names use for printing
    // the result will be "x*y + x^2"
    println(poly1.toString(Array("x", "y")))
    // the result will be "y*x + y^2"
    println(poly1.toString(Array("y", "x")))
  }

  @Test
  def test2: Unit = {
    // "x" is the first variable "y" is the second
    val ring = MultivariateRing(Z, Array("x", "y"))
    // parse polynomial
    val poly = ring("x^2 + x*y")
    // the result will be "x*y + x^2"
    println(ring show poly)
  }

  @Test
  def test3: Unit = {

    val ring = UnivariateRing(Z, "x")
    // parse polynomials
    val (p1, p2, p3) = ring("x", "x^2", "x^3")

    // this WILL modify poly1
    p1.add(p2)
    // this will NOT modify poly2
    p2.copy().add(p3)
    // this will NOT modify poly2
    ring.add(p2, p3)
  }

  @Test
  def test4: Unit = {
    import syntax._

    val ring = Zp(2)
    val a = ring("2/3")
    val b = ring("3/4")

    println(a + b)
  }

  @Test
  def test5: Unit = {
    import syntax._

    implicit val ring = UnivariateRing(Zp(3), "x")
    val (a, b) = ring("1 + 2*x^2", "1 - x")

    // add two elements
    val add = a + b
    // subtract two elements
    val sub = a - b
    // multiply two elements
    val mul = a * b
    // divide two elements
    val div = a / b
    // divide with remainder
    val divRem = a /% b
    // add one
    val inc = a ++
    // subtract one
    val dec = a --
    // negate element
    val neg = -a


  }

  @Test
  def test6: Unit = {
    import syntax._

    val a: Integer = 10
    val b: Integer = 11

    // compiles to a.add(b) (integer addition)
    assert(a + b === 21)


    implicit val ring = Zp(13)
    // compiles to ring.add(a, b) (addition mod 13)
    assert(a + b === 8)
  }

  @Test
  def test7: Unit = {
    import syntax._

    implicit val ring = UnivariateRingZp64(17, "x")
    // some random divider
    val divider = ring.randomElement()
    // some random dividend
    val dividend = 1 + 2 * divider + 3 * divider.pow(2)

    // quotient and remainder using built-in methods
    val (divPlain, remPlain) = dividend /% divider

    // precomputed Newton inverses, need to calculate it only once
    implicit val invMod: PrecomputedInverse[UnivariatePolynomialZp64] = divider.precomputedInverses
    // quotient and remainder computed using fast
    // algorithm with precomputed Newton inverses
    val (divFast, remFast) = dividend /%% divider

    // result is the same
    assert((divPlain, remPlain) == (divFast, remFast))
  }

  @Test
  def test8: Unit = {
    import rings.poly.univar.UnivariateGCD._
    import syntax._

    // Polynomials over field
    val ringZp = UnivariateRingZp64(17, "x")
    val a = ringZp("1 + 3*x + 2*x^2")
    val b = ringZp("1 - x^2")
    // Euclid and Half-GCD algorithms for polynomials over field
    assert(EuclidGCD(a, b) == HalfGCD(a, b))
    // Extended Euclidean algorithm
    val (gcd, s, t) = ExtendedEuclidGCD(a, b) match {
      case Array(gcd, s, t) => (gcd, s, t)
    }
    assert(a * s + b * t == gcd)
    // Extended Half-GCD algorithm
    val (gcd1, s1, t1) = ExtendedHalfGCD(a, b) match {
      case Array(gcd, s, t) => (gcd, s, t)
    }
    assert((gcd1, s1, t1) == (gcd, s, t))


    // Polynomials over Z
    val ringZ = UnivariateRing(Z, "x")
    val aZ = ringZ("1 + 3*x + 2*x^2")
    val bZ = ringZ("1 - x^2")
    // GCD for polynomials over Z
    assert(ModularGCD(aZ, bZ) == ringZ("1 + x"))

    // Bivariate polynomials represented as Z[y][x]
    val ringXY = UnivariateRing(UnivariateRing(Z, "y"), "x")
    val aXY = ringXY("(1 + y) + (1 + y^2)*x + (y - y^2)*x^2")
    val bXY = ringXY("(3 + y) + (3 + 2*y + y^2)*x + (3*y - y^2)*x^2")
    // Subresultant sequence
    val subResultants = SubresultantRemainders(aXY, bXY)
    // The GCD
    val gcdXY = subResultants.gcd.primitivePart
    assert(aXY % gcdXY === 0 && bXY % gcdXY === 0)
  }


  @Test
  def test9: Unit = {
    import rings.poly.PolynomialMethods._
    import rings.poly.univar.UnivariateSquareFreeFactorization._
    import rings.scaladsl.syntax._

    // ring GF(13^5)[x] (coefficient domain is finite field)
    val ringF = UnivariateRing(GF(13, 5, "z"), "x")
    // some random polynomial composed from some factors
    val polyF = ringF.randomElement() * ringF.randomElement() * ringF.randomElement().pow(10)
    // perform square-free factorization
    println(ringF show SquareFreeFactorization(polyF))
    // perform complete factorization
    println(ringF show Factor(polyF))


    // ring Q[x]
    val ringQ = UnivariateRing(Q, "x")
    // some random polynomial composed from some factors
    val polyQ = ringQ.randomElement() * ringQ.randomElement() * ringQ.randomElement().pow(10)
    // perform square-free factorization
    println(ringQ show SquareFreeFactorization(polyQ))
    // perform complete factorization
    println(ringQ show Factor(polyQ))
  }


  @Test
  def test10: Unit = {
    import rings.poly.univar.IrreduciblePolynomials._
    val random = new Random()

    // random irreducible polynomial in Z/2[x] of degree 10
    val poly1: UnivariatePolynomialZp64 = randomIrreduciblePolynomial(2, 10, random)
    assert(poly1.degree() == 10)
    assert(irreducibleQ(poly1))

    // random irreducible polynomial in Z/2[x] of degree 10
    val poly2: UnivariatePolynomial[Integer] = randomIrreduciblePolynomial(Zp(2).theRing, 10, random)
    assert(poly2.degree() == 10)
    assert(irreducibleQ(poly2))

    // random irreducible polynomial in GF(11^15)[x] of degree 10
    val poly3: UnivariatePolynomial[UnivariatePolynomialZp64] = randomIrreduciblePolynomial(GF(11, 15).theRing, 10, random)
    assert(poly3.degree() == 10)
    assert(irreducibleQ(poly3))

    // random irreducible polynomial in Z[x] of degree 10
    val poly4: UnivariatePolynomial[Integer] = randomIrreduciblePolynomialOverZ(10, random)
    assert(poly4.degree() == 10)
    assert(irreducibleQ(poly4))
  }

  @Test
  def test11: Unit = {
    import rings.poly.univar.UnivariateInterpolation._

    // points
    val points = Array(1L, 2L, 3L, 12L)
    // values
    val values = Array(3L, 2L, 1L, 6L)

    val result = new InterpolationZp64(Zp64(17))
      .update(points, values)
      .getInterpolatingPolynomial

    assert(points.zipWithIndex.forall { case (point, i) => result.evaluate(point) == values(i) })

  }


  @Test
  def test12: Unit = {

    /*  Lagrange interpolation formula */
    def lagrange[Poly <: IUnivariatePolynomial[Poly], E](points: Seq[E], values: Seq[E])(implicit ring: IUnivariateRing[Poly, E]) = {
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

    import rings.poly.univar.UnivariateInterpolation._
    import syntax._

    // coefficient ring GF(13, 5)
    implicit val cfRing = GF(13, 5, "z")
    val z = cfRing("z")
    // some points
    val points = Array(1 + z, 2 + z, 3 + z, 12 + z)
    // some values
    val values = Array(3 + z, 2 + z, 1 + z, 6 + z)

    // interpolate with Newton iterations
    val withNewton = new Interpolation(cfRing)
      .update(points, values)
      .getInterpolatingPolynomial
    // interpolate using Lagrange formula
    val withLagrange = lagrange(points, values)(UnivariateRing(cfRing, "x"))
    // results are the same
    assert(withNewton == withLagrange)
  }
}