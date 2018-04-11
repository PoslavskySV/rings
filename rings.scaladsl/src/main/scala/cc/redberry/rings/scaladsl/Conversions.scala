package cc.redberry.rings.scaladsl

import cc.redberry.rings.poly

import scala.collection.mutable.ListBuffer

/**
  *
  * @since 2.3
  */
object Conversions {

  import cc.redberry.rings.poly.multivar.MultivariateConversions

  def split[Poly <: AMultivariatePolynomial[_, Poly]](p: Poly, vars: Int*)
  : MultivariatePolynomial[Poly] =
    MultivariateConversions.split[Poly](p, vars: _*)

  def split[Poly <: AMultivariatePolynomial[_, Poly], Cf](r: IPolynomialRing[Poly, Cf], vars: Int*)
  : MultivariateRing[Poly] = {
    val raw = MultivariateConversions.split[Poly](r, vars: _*)
    val variables = r.variables.toSeq
    val cfVars = vars map variables
    val maVars = variables diff cfVars
    MultivariateRing(PolynomialRing(raw.factory().ring.getZero).setVariableNames(cfVars.toArray), maVars.toArray)
  }

  def asUnivariate[Poly <: AMultivariatePolynomial[_, Poly]](p: Poly, variable: Int)
  : UnivariatePolynomial[Poly] =
    MultivariateConversions.asUnivariate[Poly](p, variable)

  def asUnivariate[Poly <: AMultivariatePolynomial[_, Poly], Cf](r: IPolynomialRing[Poly, Cf], variable: Int)
  : UnivariateRing[Poly] = {
    val raw: poly.IPolynomialRing[UnivariatePolynomial[Poly]] = MultivariateConversions.asUnivariate[Poly](r, variable)
    val cfVars = ListBuffer[String](r.variables: _*)
    cfVars.remove(variable)
    UnivariateRing(PolynomialRing(raw.factory().ring.getZero).setVariableNames(cfVars.toArray), r.variables(variable))
  }

  def fromUnivariate[Poly <: AMultivariatePolynomial[_, Poly]](p: UnivariatePolynomial[Poly], variable: Int)
  : Poly = MultivariateConversions.fromUnivariate[Poly](p, variable)

  def fromUnivariate[Poly <: AMultivariatePolynomial[_, Poly], Cf](r: UnivariateRing[Poly], variable: Int)
  : IPolynomialRing[Poly, Cf] = {
    val raw: poly.IPolynomialRing[Poly] = MultivariateConversions.fromUnivariate[Poly](r, variable)
    val vars = r.cfRing.asInstanceOf[IPolynomialRing[_, _]].variables.toSeq
    PolynomialRing(raw.factory()).setVariableNames((vars.take(variable) ++ Seq(r.variable) ++ vars.drop(variable)).toArray)
  }
}
