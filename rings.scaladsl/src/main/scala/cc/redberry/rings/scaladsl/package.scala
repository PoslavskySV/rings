package cc.redberry.rings


/**
  * @since 1.0
  */
package object scaladsl extends Predef {

  type Ordering = java.util.Comparator[DegreeVector[_]]
  type IntZ = bigint.BigInteger

  type IPolynomial[P <: poly.IPolynomial[P]]
  = poly.IPolynomial[P]
  type IUnivariatePolynomial[P <: poly.univar.IUnivariatePolynomial[P]]
  = poly.univar.IUnivariatePolynomial[P]

  type UnivariatePolynomial[E] = poly.univar.UnivariatePolynomial[E]
  type UnivariatePolynomialZp64 = poly.univar.UnivariatePolynomialZp64
  type PrecomputedInverse[Poly <: IUnivariatePolynomial[Poly]] = poly.univar.UnivariateDivision.InverseModMonomial[Poly]

  type DegreeVector[E <: poly.multivar.DegreeVector[E]] = poly.multivar.DegreeVector[E]
  type MonomialZp64 = poly.multivar.MonomialZp64
  type Monomial[E] = poly.multivar.Monomial[E]

  type AMultivariatePolynomial[
  T <: poly.multivar.DegreeVector[T],
  P <: poly.multivar.AMultivariatePolynomial[T, P]]
  = poly.multivar.AMultivariatePolynomial[T, P]
  type MultivariatePolynomial[E] = poly.multivar.MultivariatePolynomial[E]
  type MultivariatePolynomialZp64 = poly.multivar.MultivariatePolynomialZp64

  type FactorDecomposition[P <: IPolynomial[P]] = poly.FactorDecomposition[P]

  private[scaladsl] trait LowPrioritySyntax
    extends PolynomialSetSyntax
      with PolynomialCfSyntax
      with UnivariateSyntax
      with UnivariateCfSyntax
      with MultivariateSyntax
      with MultivariateCfSyntax
      with IntegerSyntax

  object syntax extends LowPrioritySyntax {
    implicit def cfOps[E, Poly <: IPolynomial[Poly]](self: E)(implicit pRing: IPolynomialRing[Poly, E])
    = new CfOps[E, Poly](self)(pRing)

    implicit def ringOps[E](lhs: E)(implicit ringSupport: RingSupport[E]): RingOps[E] = new RingOps[E](lhs)(ringSupport.ringEv(lhs))
  }
}
