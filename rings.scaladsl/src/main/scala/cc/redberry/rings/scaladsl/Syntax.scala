package cc.redberry.rings.scaladsl

import scala.language.{implicitConversions, postfixOps}

trait RingSyntax {
  implicit def ringOps[E](lhs: E)(implicit ringSupport: RingSupport[E]): RingOps[E] = new RingOps[E](lhs)(ringSupport.ringEv(lhs))
}

trait PolynomialCfSyntax {
  implicit def polynomialOps[Poly <: IPolynomial[Poly], E](lhs: Poly)(implicit pRing: IPolynomialRing[Poly, E]): PolynomialCfOps[Poly, E] = new PolynomialCfOps[Poly, E](lhs)(pRing)
}

trait UnivariateSyntax {
  implicit def univariateOps[Poly <: IUnivariatePolynomial[Poly]](lhs: Poly)(implicit ringSupport: RingSupport[Poly]): UnivariateOps[Poly] = new UnivariateOps[Poly](lhs)(ringSupport.ringEv(lhs))
}

trait UnivariateCfSyntax {
  implicit def univariateCfOps[Poly <: IUnivariatePolynomial[Poly], E](lhs: Poly)(implicit ring: IUnivariateRing[Poly, E])
  : UnivariateCfOps[Poly, E] = new UnivariateCfOps[Poly, E](lhs)(ring)
}

trait MultivariateSyntax {
  implicit def multivariateOps[Poly <: AMultivariatePolynomial[_, Poly]]
  (lhs: Poly)(implicit ringSupport: RingSupport[Poly]): MultivariateOps[Poly]
  = new MultivariateOps[Poly](lhs)(ringSupport.ringEv(lhs))

  implicit def multivariateTermOps[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly]]
  (lhs: Poly)(implicit ringSupport: RingSupport[Poly]): MultivariateTermOps[Term, Poly]
  = new MultivariateTermOps[Term, Poly](lhs)(ringSupport.ringEv(lhs))

  implicit def multivariateTermOpsZp64
  (lhs: MultivariatePolynomialZp64)(implicit ringSupport: RingSupport[MultivariatePolynomialZp64])
  : MultivariateTermOps[MonomialZp64, MultivariatePolynomialZp64]
  = new MultivariateTermOps[MonomialZp64, MultivariatePolynomialZp64](lhs)(ringSupport.ringEv(lhs))

  implicit def multivariateTermOpsE[E]
  (lhs: MultivariatePolynomial[E])(implicit ringSupport: RingSupport[MultivariatePolynomial[E]])
  : MultivariateTermOps[Monomial[E], MultivariatePolynomial[E]]
  = new MultivariateTermOps[Monomial[E], MultivariatePolynomial[E]](lhs)(ringSupport.ringEv(lhs))
}

trait MultivariateCfSyntax {
  implicit def multivariateCfOps[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (lhs: Poly)(implicit ring: IMultivariateRing[Term, Poly, E]): MultivariateCfOps[Term, Poly, E]
  = new MultivariateCfOps[Term, Poly, E](lhs)(ring)

  implicit def multivariateCfOpsZp64
  (lhs: MultivariatePolynomialZp64)(implicit ring: IMultivariateRing[MonomialZp64, MultivariatePolynomialZp64, Long])
  : MultivariateCfOps[MonomialZp64, MultivariatePolynomialZp64, Long]
  = new MultivariateCfOps[MonomialZp64, MultivariatePolynomialZp64, Long](lhs)(ring)

  implicit def multivariateCfOpsE[E]
  (lhs: MultivariatePolynomial[E])(implicit ring: IMultivariateRing[Monomial[E], MultivariatePolynomial[E], E])
  : MultivariateCfOps[Monomial[E], MultivariatePolynomial[E], E]
  = new MultivariateCfOps[Monomial[E], MultivariatePolynomial[E], E](lhs)(ring)

  def multivariateImplicits[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (implicit ring: IMultivariateRing[Term, Poly, E])
  : Poly => MultivariateCfOps[Term, Poly, E] = p => new MultivariateCfOps[ring.MonomialType, ring.ElementType, ring.CoefficientType](p)(ring)
}

trait CfSyntax {
  implicit def cfOps[E, Poly <: IPolynomial[Poly]](self: E)(implicit pRing: IPolynomialRing[Poly, E])
  = new CfOps[E, Poly](self)(pRing)
}

trait IntegerSyntax {
  implicit def integerOps[E](self: Int)(implicit ring: Ring[E]) = new IntegerOps[E](self)(ring)
}