package cc.redberry.rings

import java.util.Comparator

import cc.redberry.rings.bigint.BigInteger

/**
  *
  * @since 1.0
  */
package object scaladsl extends Predef {

  type Ordering = Comparator[DegreeVector[_]]
  type Integer = BigInteger

  type IPolynomial[P <: cc.redberry.rings.poly.IPolynomial[P]]
  = cc.redberry.rings.poly.IPolynomial[P]
  type IUnivariatePolynomial[P <: cc.redberry.rings.poly.univar.IUnivariatePolynomial[P]]
  = cc.redberry.rings.poly.univar.IUnivariatePolynomial[P]

  type UnivariatePolynomial[E] = cc.redberry.rings.poly.univar.UnivariatePolynomial[E]
  type UnivariatePolynomialZp64 = cc.redberry.rings.poly.univar.UnivariatePolynomialZp64

  type DegreeVector[E <: cc.redberry.rings.poly.multivar.DegreeVector[E]] = cc.redberry.rings.poly.multivar.DegreeVector[E]
  type MonomialZp64 = cc.redberry.rings.poly.multivar.MonomialZp64
  type Monomial[E] = cc.redberry.rings.poly.multivar.Monomial[E]

  type AMultivariatePolynomial[
  T <: cc.redberry.rings.poly.multivar.DegreeVector[T],
  P <: cc.redberry.rings.poly.multivar.AMultivariatePolynomial[T, P]]
  = cc.redberry.rings.poly.multivar.AMultivariatePolynomial[T, P]
  type MultivariatePolynomial[E] = cc.redberry.rings.poly.multivar.MultivariatePolynomial[E]
  type MultivariatePolynomialZp64 = cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64

  type FactorDecomposition[P <: IPolynomial[P]] = cc.redberry.rings.poly.FactorDecomposition[P]

  private[scaladsl] trait LowPrioritySyntax
    extends PolynomialSetSyntax
      with PolynomialCfSyntax
      with UnivariateSyntax
      with UnivariateCfSyntax
      with MultivariateSyntax
      with MultivariateCfSyntax
      with IntegerSyntax

  object syntax extends LowPrioritySyntax {
    implicit def cfOps[E, Poly <: IPolynomial[Poly]](self: E)(implicit pRing: PolynomialRing[Poly, E])
    = new CfOps[E, Poly](self)(pRing)

    implicit def ringOps[E](lhs: E)(implicit ringSupport: RingSupport[E]): RingOps[E] = new RingOps[E](lhs)(ringSupport.ringEv(lhs))
  }
}
