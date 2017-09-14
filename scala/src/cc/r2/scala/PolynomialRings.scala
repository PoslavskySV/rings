package cc.r2.scala

import java.util.Comparator

import cc.r2.core.number.BigInteger
import cc.r2.core.poly._
import cc.r2.core.poly.multivar.{DegreeVector, MultivariatePolynomial, MultivariatePolynomialZp64, lMultivariatePolynomialZp}
import cc.r2.core.poly.univar.{UnivariatePolynomial, UnivariatePolynomialZp64}

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
object PolynomialRings {
  abstract class PolynomialRing[Poly <: IPolynomial[Poly]]
  (val pDomain: PolynomialDomain[Poly], val variables: Array[String]) {
    final def parse(string: String): Poly = pDomain.parse(string, variables)

    final def show(polynomial: Poly): String = polynomial.toString(variables: _*)

    final def index(variable: String) = variables.indexOf(variable)
  }

  implicit def domainMethods[Poly <: IPolynomial[Poly]](ring: PolynomialRing[Poly]) = ring.pDomain

  final case class GaloisField64
  (override val pDomain: FiniteField[UnivariatePolynomialZp64], variable: String)
    extends PolynomialRing[UnivariatePolynomialZp64](pDomain, Array(variable))

  object GaloisField64 {
    def apply(modulus: Long, exponent: Int, variable: String) = new GaloisField64(Domains.GF(modulus, exponent), variable)

    def apply(irreducible: UnivariatePolynomialZp64, variable: String) = new GaloisField64(Domains.GF(irreducible), variable)
  }

  final case class GaloisField[E]
  (override val pDomain: FiniteField[UnivariatePolynomial[E]], variable: String)
    extends PolynomialRing[UnivariatePolynomial[E]](pDomain, Array(variable))

  object GaloisField {
    def apply(modulus: BigInteger, exponent: Int, variable: String) = new GaloisField[BigInteger](Domains.GF(modulus, exponent), variable)

    def apply[E](irreducible: UnivariatePolynomial[E], variable: String) = new GaloisField[E](Domains.GF(irreducible), variable)
  }

  final case class UnivariateRingZp64
  (modulus: IntegersZp64, variable: String)
    extends PolynomialRing[UnivariatePolynomialZp64](Domains.PolynomialsZp(modulus), Array(variable))

  object UnivariateRingZp64 {
    def apply(modulus: Long, variable: String) = new UnivariateRingZp64(new IntegersZp64(modulus), variable)
  }

  final case class UnivariateRing[Element]
  (domain: Domain[Element], variable: String)
    extends PolynomialRing[UnivariatePolynomial[Element]](Domains.Polynomials(domain), Array(variable))

  private type MonomialOrder = Comparator[DegreeVector[_]]

  final case class MultivariateRingZp64
  (modulus: IntegersZp64, override val variables: Array[String], ordering: MonomialOrder)
    extends PolynomialRing[MultivariatePolynomialZp64](Domains.MultivariatePolynomialsZp(variables.length, modulus), variables) {
    def this(modulus: Long, variables: Array[String], ordering: MonomialOrder) = this(new IntegersZp64(modulus), variables, ordering)
  }

  final case class MultivariateRing[ElementType]
  (domain: Domain[ElementType], override val variables: Array[String], ordering: MonomialOrder)
    extends PolynomialRing[MultivariatePolynomial[ElementType]](Domains.MultivariatePolynomials(variables.length, domain), variables)
}
