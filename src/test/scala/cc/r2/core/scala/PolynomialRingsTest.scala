package cc.r2.core.scala

import cc.r2.core.poly.PolynomialMethods
import cc.r2.core.poly.multivar.MonomialOrder
import cc.r2.core.poly.univar.{IrreduciblePolynomials, UnivariatePolynomialZ64}
import cc.r2.core.scala.PolynomialRings.{GF, MultivariateRing, UnivariateRing, UnivariateRingZp64}
import org.apache.commons.math3.random.Well1024a
import org.junit.Test

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
class PolynomialRingsTest {
  @Test
  def simpleUnivariateRing = {
    import implicits._

    implicit val ring = UnivariateRingZp64(17, "z")
    println(s"Ring $ring")

    // parse univariate polynomial from string
    val poly1 = parse("1 + 2*z^2 - z^3 + z^5 - z^17")

    // create univariate polynomial with factory method
    val poly2 = UnivariatePolynomialZ64.create(0, 1, 2, 3, 4).modulus(ring.coefficientDomain)

    // another way to create poly
    val x = parse("z")
    val poly3 = 1 + (x ** 3) + 3 * x ** 9 + 9 * x ** 27

    // mix polynomials in some way
    val poly4 = (poly2 - poly1) ** 2 + 2L * (poly2 - poly1) * (poly3 - poly1) + (poly3 - poly1) ** 2

    // calculate GCD
    val gcd = PolynomialMethods.PolynomialGCD(poly1 * poly3, poly2 * poly4)
    println(s"gcd : ${ring show gcd}")

    // factor some complicated poly
    val factors = PolynomialMethods.factor(poly1 * poly2 ** 2 * poly3 ** 3 * poly4 ** 4)
    println(s"factors : ${ring show factors}")
  }

  @Test
  def compicatedGaloisField = {
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
    println(PolynomialMethods.factor(rndPoly))
  }

  @Test
  def complicatedMultivariateRing = {
    import implicits._

    val coefficientRing = GF(UnivariateRingZp64(41, "x")("1 + x + x^2"), "z")

    implicit val ring = MultivariateRing(coefficientRing, Array("a", "b", "c"), MonomialOrder.LEX)

    println(s"Coefficient ring: $coefficientRing")

    // parse univariate polynomial from string
    val poly1 = parse("(1 + x + x^2)*a*b*c + a^2 + (x)*c^3 + (13 + 12*x + 4*x^6)*b^2*c^3 + (2 + x)")
    val poly2 = parse("(1 - x - x^2)*a*b*c - a^2 - (x)*c^3 + (13 - 12*x + 4*x^6)*b^2*c^3 + (2 + x^87)")
    val poly = poly1 * poly2 ** 2 * (poly1 + poly2) ** 3

    println(PolynomialMethods.factor(poly))
  }
}
