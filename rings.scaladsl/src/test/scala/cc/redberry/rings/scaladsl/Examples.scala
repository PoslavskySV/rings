package cc.redberry.rings.scaladsl

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.multivar._
import cc.redberry.rings.poly.univar._
import cc.redberry.rings.primes.BigPrimes
import org.junit.Test

import scala.runtime.ScalaRunTime

/**
  *
  * @since 1.0
  */
class Examples {
  @Test
  def factorizationInZ17: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    // Z/17[x]
    implicit val ring = UnivariateRingZp64(17, "x")
    println(s"Ring $ring")

    // parse univariate polynomial from string
    val poly1 = ring("1 + 2*x + 3*x^2")
    val poly2 = ring("4 + 5*x^5 + 6*x^6")
    val poly3 = ring("7 + 8*x^9 + 9*x^10 + 10*x^11 + x^99")
    val poly = poly1 * poly2 * poly3

    // factorize poly
    val factors = Factor(poly)
    // there will be exactly 9 factors
    assert(9 == factors.size())

    println(s"factorization : ${ring show factors}")
  }


  @Test
  def factorizationInZ1237940039285380274899124357: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    // prime = 1237940039285380274899124357
    val prime = BigPrimes.nextPrime(BigInteger.valueOf(1).shiftLeft(100))
    val coefficientRing = Zp(prime)
    // Z/1237940039285380274899124357[x]
    implicit val ring = UnivariateRing(coefficientRing, "x")

    println(s"Ring $ring")

    // parse univariate polynomial from string
    val poly1 = ring("1 + 2*x + 3*x^2")
    val poly2 = ring("4 + 5*x^5 + 6*x^6")
    val poly3 = ring("7 + 8*x^9 + 9*x^10 + 10*x^11 + x^99")
    val poly = poly1 * poly2 * poly3

    // factorize poly
    val factors = Factor(poly)
    // there will be exactly 9 factors
    assert(9 == factors.size())

    println(s"factorization : ${ring show factors}")
  }

  @Test
  def factorizationInGFxxxx: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    val cfRing = GF(7, 10, "z")
    val z = cfRing("z")
    // Z/17[x]
    implicit val ring = MultivariateRing(cfRing, Array("a", "b", "c"))
    println(s"Ring $ring")

    val (a, b, c) = ring("a", "b", "c")

    val xxxx = a**2 + 2 + b + c**123123

    // parse univariate polynomial from string
    val poly1: MultivariatePolynomial[UnivariatePolynomialZp64] = ring("(1 + z)*a*b^2 + c^5")
    val poly2 = ring("(1 + 5*z + z^67)*a^6*b^2*c^8 + (z)*a*b*c^5 + 1")
    val poly3 = ring("2131 + a^7 + b^6 + c^7")
    val poly = poly1 * poly2 * poly3

    val lc: UnivariatePolynomialZp64 = poly1.lc
    val xxx = lc ** 2 + 1

    // factorize poly
    val factors = Factor(poly)
    //    println(factors)

    println(poly1)
    println(poly2)
    println(poly3)
    println(s"factorization : ${ring show factors}")
  }

  @Test
  def factorizationInZ2xyz: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    // Z/2[x,y,z]
    implicit val ring = MultivariateRingZp64(2, Array("x", "y", "z"), MonomialOrder.LEX)

    val (x, y, z) = ring("x", "y", "z")
    val poly1 = 1 + (1 + x + y + z) ** 3
    val poly2 = x + (1 + x + y) ** 3
    val poly3 = y + (1 + x + z) ** 3
    val poly4 = z + (1 + y + z) ** 3
    val poly = poly1 * poly2 * poly3 * poly4

    // factorize polynomial
    val factors = Factor(poly)
    // there will be exactly 5 factors
    assert(5 == factors.size())

    println(s"factorization : ${ring show factors}")
  }

  @Test
  def factorizationInZxyz: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    // Z/2[x,y,z]
    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"), MonomialOrder.LEX)
    import ring.{`x`, `y`, `z`}

    val poly1 = 1 + ring("1 + x + y + z") ** 3
    val poly2 = `x` + ring("1 + x + y") ** 3
    val poly3 = `y` + ring("1 + x + z") ** 3
    val poly4 = `z` + ring("1 + y + z") ** 3
    val poly = poly1 * poly2 * poly3 * poly4

    // factorize polynomial
    val factors = Factor(poly)
    // there will be exactly 5 factors
    assert(5 == factors.size())

    println(s"factorization : ${ring show factors}")
  }

  @Test
  def gcdInZxyz: Unit = {

    import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.scaladsl.Rings._
    import cc.redberry.rings.scaladsl.implicits._

    // Z/2[x,y,z]
    implicit val ring = MultivariateRing(Q, Array("x1", "x2", "x3", "x4", "x5"), MonomialOrder.LEX)

    // some poly
    val a = ring("(2/3)*x1^6*x2 + (1/7)*x1*x2*x4 + 2")
    val b = ring("1 + (2/7)*x2^16 + x3^34 + x5^16 + (2/3)*x1*x2*x3*x4*x5")
    val gcd = a + b

    val actualGCD = PolynomialGCD(
      ((a - b) ** 2 + b ** 4 + a ** 2) * gcd,
      (b - a + a ** 2) * gcd)

    assert(gcd.monic() == actualGCD.monic())

    println(s"GCD : ${ring show gcd}")

    println(MultivariateGCD.EEZGCD[ring.MonomialType, ring.ElementType](a, b))
  }


  @Test
  def rings: Unit = {

    import cc.redberry.rings.scaladsl.Rings._

    // Ring Z[x]
    UnivariateRing(Z, "x")
    // Ring Z[x, y, z]
    MultivariateRing(Z, Array("x", "y", "z"))
    // Ring Q[a, b, c]
    MultivariateRing(Q, Array("a", "b", "c"))

    // Ring Z/3[x] (64 indicates that machine numbers are used in the basis)
    UnivariateRingZp64(3, "x")
    // Ring Z/3[x, y, z]
    MultivariateRingZp64(3, Array("x", "y", "z"))
    // Ring Z/1267650600228229401496703205653[x, y, z]
    MultivariateRing(Zp((BigInt(2) ^ 107) - 1), Array("x", "y", "z"))

    // Galois field with cardinality 7^10
    GF(7, 10, "x")

    ScalaRunTime.stringOf()
    Rationals(UnivariateRing(Z, "x"))
  }

}

