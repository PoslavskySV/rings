package cc.redberry.rings.scaladsl

import cc.redberry.rings.poly.multivar.MultivariatePolynomial
import org.junit.Test

import scala.language.postfixOps

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
}