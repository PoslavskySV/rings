package cc.r2.core.scala

import java.util.Comparator

import cc.r2.core.number.{BigInteger, Rational}
import cc.r2.core.poly.multivar.{DegreeVector, MonomialOrder, MultivariatePolynomial, MultivariatePolynomialZp64}
import cc.r2.core.poly.univar.{UnivariatePolynomial, UnivariatePolynomialZp64}
import cc.r2.core.poly.{ElementParser, _}

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
object PolynomialRings {

  /**
    * Simple wrapper around [[Domain]] used to unify PolynomialRing and Domain
    *
    * @param domain the [[Domain]]
    **/
  sealed class Ring[E](val domain: Domain[E]) extends ToStringSupport[E] with ElementParser[E] {
    /**
      * @inheritdoc
      **/
    override def toString(element: E): String = element.toString

    /**
      * @inheritdoc
      **/
    override def parse(string: String): E = domain parse string

    /**
      * @inheritdoc
      **/
    override def toString: String = domain.toString

    private val eClass: Class[_ <: E] = domain.getZero.getClass

    /**
      * Reflection: determines whether is element of this
      */
    final def isElement(e: Any): Boolean = eClass.isAssignableFrom(e.getClass)

    /**
      * Casts e to element of this
      */
    final def element(e: Any): E = e.asInstanceOf[E]
  }

  /**
    * Delegate [[Domain]] methods for [[Ring]]
    */
  implicit def domainMethods[E](ring: Ring[E]): Domain[E] = ring.domain

  /**
    * Implicitly convert [[Domain]] to [[Ring]]
    */
  implicit def domainAsRing[E](domain: Domain[E]) = new Ring[E](domain)

  object Rings {
    val Z: Ring[BigInteger] = Domains.Z

    val Q: Ring[Rational] = Domains.Q

    def Zp(modulus: Long): Ring[BigInteger] = Domains.Zp(BigInteger.valueOf(modulus))

    def Zp(modulus: BigInteger): Ring[BigInteger] = Domains.Zp(modulus)
  }

  /**
    * Base class for polynomial rings
    *
    * @param domain    the [[PolynomialDomain]]
    * @param variables polynomial variables
    */
  sealed abstract class PolynomialRing[Poly <: IPolynomial[Poly]]
  (override val domain: PolynomialDomain[Poly],
   val variables: Array[String]) extends Ring[Poly](domain) {

    /**
      * Parse polynomial
      */
    override def parse(string: String): Poly = domain.parse(string, variables)

    /**
      * Parse polynomial
      */
    final def apply(string: String) = parse(string)

    /**
      * String representation of polynomial from this ring
      */
    def show(obj: WithVariables): String = obj match {
      case fct: FactorDecomposition[Poly] => fct.toString(this)
      case _ => obj.toString(variables)
    }

    /**
      * String representation of a seq of polynomials
      */
    final def show(list: Traversable[_]): String
    = list.map {
      case wv: WithVariables => show(wv)
      case p => p.toString
    }.toString

    /**
      * Index of variable with specified string representation
      */
    final def index(variable: String) = variables.indexOf(variable)

    /**
      * @inheritdoc
      **/
    override final def toString(element: Poly): String = show(element)

    /**
      * String representation of this ring
      */
    override def toString: String = domain.toString(variables)
  }

  /**
    * Delegate [[PolynomialDomain]] methods for [[PolynomialRing]]
    */
  implicit def domainMethods[Poly <: IPolynomial[Poly]](ring: PolynomialRing[Poly]): PolynomialDomain[Poly] = ring.domain

  /**
    * Galois field with prime base in a range of `(0, 2^63)`
    *
    * @param domain   the [[FiniteField]]
    * @param variable the variable of univariate polynomials representing this Galois field
    */
  final case class GaloisField64
  (override val domain: FiniteField[UnivariatePolynomialZp64], variable: String)
    extends PolynomialRing[UnivariatePolynomialZp64](domain, Array(variable))

  /**
    * Delegate [[FiniteField]] methods fo [[GaloisField64]]
    */
  implicit def domainMethods(ring: GaloisField64): FiniteField[UnivariatePolynomialZp64] = ring.domain

  /**
    * Galois field with arbitrary prime base
    *
    * @param domain   the [[FiniteField]]
    * @param variable the variable of univariate polynomials representing this Galois field
    */
  final case class GaloisField[E]
  (override val domain: FiniteField[UnivariatePolynomial[E]], variable: String,
   private val ringOption: Option[UnivariateRing[E]] = None)
    extends PolynomialRing[UnivariatePolynomial[E]](domain, Array(variable)) {
    /**
      * @inheritdoc
      */
    override def parse(string: String): UnivariatePolynomial[E] = ringOption match {
      case Some(ring) => domain.valueOf(domain.getOne.parsePoly(string, ring.coefficientDomain, variable))
      case None => super.parse(string)
    }

    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = ringOption match {
      case Some(ring) =>
        obj match {
          case poly: UnivariatePolynomial[E] => poly.toString(ring.coefficientDomain, variable)
          case _ => super.show(obj)
        }
      case None => super.show(obj)
    }

    /**
      * @inheritdoc
      */
    override val toString: String = ringOption match {
      case Some(ring) => domain.toString(ring.coefficientDomain.toString, this, variables)
      case None => super.toString
    }
  }

  /**
    * Delegate [[FiniteField]] methods fo [[GaloisField64]]
    */
  implicit def domainMethods[E](ring: GaloisField[E]): FiniteField[UnivariatePolynomial[E]] = ring.domain

  object GF {
    /**
      * Create Galois field with cardinality `modulus ^ exponent` represented by univariate polynomials with
      * specified variable
      */
    def apply(modulus: Long, exponent: Int, variable: String)
    = GaloisField64(Domains.GF(modulus, exponent), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply(irreducible: UnivariatePolynomialZp64, variable: String)
    = GaloisField64(Domains.GF(irreducible), variable)

    /**
      * Create Galois field with cardinality `modulus ^ exponent` represented by univariate polynomials with
      * specified variable
      */
    def apply(modulus: BigInteger, exponent: Int, variable: String)
    = GaloisField[BigInteger](Domains.GF(modulus, exponent), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply[E](irreducible: UnivariatePolynomial[E], variable: String)
    = GaloisField[E](Domains.GF(irreducible), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply[E](irreducible: UnivariatePolynomial[E], univariateRing: UnivariateRing[E], variable: String)
    = GaloisField[E](Domains.GF(irreducible), variable, Some(univariateRing))
  }

  /**
    * Ring of Zp[x] polynomials
    *
    * @param coefficientDomain coefficient domain
    * @param variable          variable
    */
  final case class UnivariateRingZp64 private(coefficientDomain: IntegersZp64, variable: String)
    extends PolynomialRing[UnivariatePolynomialZp64](Domains.PolynomialsZp(coefficientDomain), Array(variable)) {
    val modulus: Long = coefficientDomain.modulus
  }

  object UnivariateRingZp64 {
    /**
      * Zp[variable] with specified modulus
      */
    def apply(modulus: Long, variable: String): UnivariateRingZp64 = UnivariateRingZp64(new IntegersZp64(modulus), variable)

    /**
      * Zp[variable] with specified coefficient domain (Zp)
      */
    def apply(domain: IntegersZp64, variable: String) = new UnivariateRingZp64(domain, variable)
  }

  /**
    * Ring of univariate polynomials over generic domains
    *
    * @param coefficientDomain coefficient domain
    * @param variable          variable
    */
  final case class UnivariateRing[E](coefficientDomain: Ring[E], variable: String)
    extends PolynomialRing[UnivariatePolynomial[E]](Domains.Polynomials(coefficientDomain.domain), Array(variable)) {

    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = obj match {
      case poly: UnivariatePolynomial[E] => poly.toString(coefficientDomain, variable)
      case _ => super.show(obj)
    }

    /**
      * @inheritdoc
      */
    override def parse(string: String): UnivariatePolynomial[E] = domain.getOne.parsePoly(string, coefficientDomain, variable)

    /**
      * @inheritdoc
      */
    override val toString: String = domain.toString(coefficientDomain.toString, variables)
  }

  private type Ordering = Comparator[DegreeVector[_]]

  /**
    * Zp[variables] with specified modulus
    *
    * @param coefficientDomain coefficient domain
    */
  final case class MultivariateRingZp64
  (coefficientDomain: IntegersZp64, override val variables: Array[String], ordering: Ordering)
    extends PolynomialRing[MultivariatePolynomialZp64](Domains.MultivariatePolynomialsZp(variables.length, coefficientDomain), variables)

  object MultivariateRingZp64 {
    /**
      * Zp[variables] with specified modulus, variables and ordering
      */
    def apply(modulus: Long, variables: Array[String], ordering: Ordering): MultivariateRingZp64
    = MultivariateRingZp64(new IntegersZp64(modulus), variables, ordering)

    /**
      * Zp[variables] with specified modulus, variables and default ordering (LEX)
      */
    def apply(modulus: Long, variables: Array[String]): MultivariateRingZp64
    = MultivariateRingZp64(new IntegersZp64(modulus), variables, MonomialOrder.LEX)
  }

  /**
    * Ring of multivariate polynomials over generic domains
    *
    * @param coefficientDomain coefficient domain
    */
  final case class MultivariateRing[E]
  (coefficientDomain: Ring[E], override val variables: Array[String], ordering: Ordering)
    extends PolynomialRing[MultivariatePolynomial[E]](Domains.MultivariatePolynomials(variables.length, coefficientDomain), variables) {
    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = obj match {
      case poly: MultivariatePolynomial[E] => poly.toString(coefficientDomain, variables)
      case _ => super.show(obj)
    }

    /**
      * @inheritdoc
      */
    override def parse(string: String): MultivariatePolynomial[E] = domain.getOne.parsePoly(string, coefficientDomain, variables)

    /**
      * @inheritdoc
      */
    override val toString: String = domain.toString(coefficientDomain.toString, variables)
  }

  object MultivariateRing {
    /**
      * Zp[variables] with specified modulus, variables and default ordering (LEX)
      */
    def apply[E](coefficientDomain: Ring[E], variables: Array[String]): MultivariateRing[E]
    = MultivariateRing[E](coefficientDomain, variables, MonomialOrder.LEX)
  }
}
