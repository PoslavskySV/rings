package cc.r2.scala

import cc.r2.core.poly._
import cc.r2.core.poly.multivar._
import cc.r2.core.poly.univar.UnivariatePolynomialZp64

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
object Main {


  def main(args: Array[String]): Unit = {

    import PolynomialRings._
    import implicits._

    //    implicit def asPoly[
    //    Term <: DegreeVector[Term],
    //    Poly <: AMultivariatePolynomial[Term, Poly],
    //    np <: AMultivariatePolynomial[_, _]](np: np) = np.asInstanceOf[Poly]

    implicit val ring = MultivariateRing(Domains.GF(17, 3), Array("x", "y", "z"), MonomialOrder.LEX)

    val x = parse("2 * x + y + z")
    val y = parse("3 * x - y^2 - z^2")
    val z: MultivariatePolynomial[UnivariatePolynomialZp64] = x + y

    println(x.mapCoefficients(Domains.Z, p => asBigInteger(p.cc())))

    val factors = MultivariateFactorization.factor(z)
    //    println(AMultivariatePolynomial.swapVariables(x, 1, 2))
    println(factors)
    val gcd = MultivariateGCD.PolynomialGCD(x, y) * x * y

    val ddd = Domains.GF(17, 3).valueOf(2L)
    println(gcd.swapVariables("x" <> "y"))
    println(gcd)
    println("sasa")


    //    implicit val ring = new GaloisField64(3, 4, "x")
    //    val x = parse("x^2 + 1")
    //    val y = parse("x^4 + 2")
    //
    //    val zp = x + y
    //    println(zp)


    //    val const = parse("2")
    //    println(const)
    //    val inc = const.++
    //    println(inc ++)


    //    implicit val variables = Variables(Array("x", "y", "z"))
    //    val a = MultivariatePolynomial.parse("x^2 + y^2", variables.names: _*)
    //    val b = MultivariatePolynomial.parse("x^4 - y^4", variables.names: _*)
    //
    //    println(a - b)
    //

    //
    //    //    type PolyZ = UnivariatePolynomial[BigInteger]
    //
    //    val x = UnivariatePolynomial.create(Domains.Z, 0, 1)
    //    val poly = 1 + x + x ** 2 - x + x ** 3 + 2 - x
    //    println(UnivariateFactorization.factor(poly))
  }
}
