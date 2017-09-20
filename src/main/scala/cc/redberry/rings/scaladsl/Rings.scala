package cc.redberry.rings.scaladsl

import java.util.Comparator

import cc.redberry.rings
import cc.redberry.rings.bigint.{BigInteger, Rational}
import cc.redberry.rings.poly._
import cc.redberry.rings.poly.multivar._
import cc.redberry.rings.poly.univar.{IUnivariatePolynomial, UnivariatePolynomial, UnivariatePolynomialZp64}
import cc.redberry.rings.{poly, _}

/**
  *
  *
  * @since 1.0
  */
object Rings {

  /**
    * Simple wrapper around [[Ring]] used to unify PolynomialRing and Ring
    *
    * @param domain the [[Ring]]
    **/
  sealed class Ring[E](val domain: rings.Ring[E]) extends ToStringSupport[E] with ElementParser[E] {
    /**
      * Element type
      */
    final type Type = E

    /**
      * @inheritdoc
      **/
    override def toString(element: E): String = element toString

    /**
      * @inheritdoc
      **/
    override def parse(string: String): E = domain parse string

    /**
      * @inheritdoc
      **/
    override def toString: String = domain toString

    /** element class  */
    private val eClass: Class[_ <: E] = domain.getZero.getClass

    /**
      * Reflection: determines whether is element of this
      */
    final def isElement(e: Any): Boolean = eClass isAssignableFrom e.getClass

    /**
      * Casts e to element of this
      */
    final def element(e: Any): E = e.asInstanceOf[E]

    /**
      * Pretty toString
      */
    def show(obj: Any): String = obj toString
  }

  /**
    * Delegate [[Ring]] methods for [[Ring]]
    */
  implicit def domainMethods[E](ring: Ring[E]): rings.Ring[E] = ring.domain

  /**
    * Implicitly convert [[Ring]] to [[Ring]]
    */
  implicit def domainAsRing[E](domain: rings.Ring[E]) = new Ring[E](domain)

  /**
    * Ring of integers (Z)
    */
  val Z: Ring[BigInteger] = rings.Rings.Z

  /**
    * Field of rationals (Q)
    */
  val Q: Ring[Rational] = rings.Rings.Q

  /**
    * Field of integers modulo `modulus`
    *
    * @param modulus the modulus
    */
  def Zp(modulus: Long): Ring[BigInteger] = rings.Rings.Zp(BigInteger.valueOf(modulus))

  /**
    * Field of integers modulo `modulus`
    *
    * @param modulus the modulus
    */
  def Zp(modulus: BigInteger): Ring[BigInteger] = rings.Rings.Zp(modulus)

  /**
    * Base class for polynomial rings
    *
    * @param domain    the [[PolynomialRing]]
    * @param variables polynomial variables
    * @tparam E coefficient type
    */
  sealed abstract class PolynomialRing[Poly <: IPolynomial[Poly], E]
  (override val domain: poly.PolynomialRing[Poly],
   val variables: Array[String]) extends Ring[Poly](domain) {

    /**
      * Type of coefficients
      */
    final type CoefficientType = E

    /**
      * Parse polynomial
      */
    override def parse(string: String): Poly = domain.parse(string, variables)

    /**
      * Parse polynomial
      */
    final def apply(string: String) = parse(string)

    /**
      * To string
      */
    override def show(obj: Any): String = obj match {
      case wv: WithVariables => show(wv)
      case _ => obj.toString
    }

    /**
      * String representation of polynomial from this ring
      */
    def show(obj: WithVariables): String = obj match {
      case fct: FactorDecomposition[Poly] => fct.toString(this)
      case el if isElement(el) => el.toString(variables)
      case _ => obj.toString()
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
    final def index(variable: String): Int = variables.indexOf(variable)

    /**
      * @inheritdoc
      **/
    override final def toString(element: Poly): String = show(element)

    /**
      * String representation of this ring
      */
    override def toString: String = domain.toString(variables)

    /**
      * The first variable
      */
    final def `x`: Poly = domain.variable(0)

    /**
      * The second variable
      */
    final def `y`: Poly = domain.variable(1)

    /**
      * The third variable
      */
    final def `z`: Poly = domain.variable(2)

    /**
      * Shortcut for /% operation
      */
    protected[Rings] final def divRem(a: Poly, b: Poly): (Poly, Poly) = {
      val qd = domain.divideAndRemainder(a, b)
      if (qd == null)
        throw new ArithmeticException(s"not divisible with remainder: ${this show a} / ${this show b}")
      (qd(0), qd(1))
    }

    /**
      * Constant polynomial with specified value
      */
    def getConstant(value: E): Poly

    /**
      * Add coefficient ring element
      */
    def add(poly: Poly, el: E): Poly

    /**
      * Subtract coefficient ring element
      */
    def subtract(poly: Poly, el: E): Poly

    /**
      * Multiply by coefficient ring element
      */
    def multiply(poly: Poly, el: E): Poly

    /**
      * Divide by coefficient ring element
      */
    def divide(poly: Poly, el: E): Poly

    /**
      * Divide by coefficient ring element
      */
    def divideAndRemainder(poly: Poly, el: E): (Poly, Poly)

    /**
      * Value of integer in coefficient ring
      */
    def cfValue(i: Int): E
  }

  /**
    * Delegate [[PolynomialRing]] methods for [[PolynomialRing]]
    */
  implicit def domainMethods[Poly <: IPolynomial[Poly], E](ring: PolynomialRing[Poly, E]): poly.PolynomialRing[Poly] = ring.domain

  /**
    * Ring of univariate polynomials
    *
    * @param domain   the [[PolynomialRing]]
    * @param variable the variable
    */
  sealed abstract class IUnivariateRing[Poly <: IUnivariatePolynomial[Poly], E]
  (override val domain: poly.PolynomialRing[Poly],
   val variable: String) extends PolynomialRing[Poly, E](domain, Array(variable))

  /**
    * Galois field with prime base in a range of `(0, 2^63)`
    *
    * @param domain   the [[FiniteField]]
    * @param variable the variable of univariate polynomials representing this Galois field
    */
  final case class GaloisField64
  (override val domain: FiniteField[UnivariatePolynomialZp64], override val variable: String)
    extends IUnivariateRing[UnivariatePolynomialZp64, Long](domain, variable) {

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: Long): UnivariatePolynomialZp64 = domain.valueOf(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().multiply(el)

    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().divide(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: UnivariatePolynomialZp64, el: Long): (UnivariatePolynomialZp64, UnivariatePolynomialZp64)
    = divRem(poly, domain.valueOf(el))

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): Long = i.asInstanceOf[Long]
  }

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
  (override val domain: FiniteField[UnivariatePolynomial[E]], override val variable: String,
   private val ringOption: Option[UnivariateRing[E]] = None)
    extends IUnivariateRing[UnivariatePolynomial[E], E](domain, variable) {

    val polyCoefficientDomain = domain.factory().ring

    /**
      * @inheritdoc
      */
    override def parse(string: String): UnivariatePolynomial[E] = ringOption match {
      case Some(ring) => domain.valueOf(domain.factory().parsePoly(string, ring.coefficientDomain, variable))
      case None => super.parse(string)
    }

    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = ringOption match {
      case Some(ring) =>
        obj match {
          case poly: UnivariatePolynomial[E] => poly.toString(ring.coefficientDomain, variable)
          case cfx if ringOption.isDefined && ringOption.get.isElement(cfx) => ringOption.get.show(cfx)
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

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: E): UnivariatePolynomial[E] = domain.factory().createConstant(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().multiply(el)


    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().divideExact(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: UnivariatePolynomial[E], el: E): (UnivariatePolynomial[E], UnivariatePolynomial[E])
    = divRem(poly, poly.createConstant(el))

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): E = polyCoefficientDomain.valueOf(i.asInstanceOf[Long])
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
    = GaloisField64(rings.Rings.GF(modulus, exponent), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply(irreducible: UnivariatePolynomialZp64, variable: String)
    = GaloisField64(rings.Rings.GF(irreducible), variable)

    /**
      * Create Galois field with cardinality `modulus ^ exponent` represented by univariate polynomials with
      * specified variable
      */
    def apply(modulus: BigInteger, exponent: Int, variable: String)
    = GaloisField[BigInteger](rings.Rings.GF(modulus, exponent), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply[E](irreducible: UnivariatePolynomial[E], variable: String)
    = GaloisField[E](rings.Rings.GF(irreducible), variable)

    /**
      * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
      * specified variable
      */
    def apply[E](irreducible: UnivariatePolynomial[E], univariateRing: UnivariateRing[E], variable: String)
    = GaloisField[E](rings.Rings.GF(irreducible), variable, Some(univariateRing))
  }

  /**
    * Ring of Zp[x] polynomials
    *
    * @param coefficientDomain coefficient ring
    * @param variable          variable
    */
  final case class UnivariateRingZp64 private(coefficientDomain: IntegersZp64, override val variable: String)
    extends IUnivariateRing[UnivariatePolynomialZp64, Long](rings.Rings.UnivariateRingZp(coefficientDomain), variable) {
    val modulus: Long = coefficientDomain.modulus

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: Long): UnivariatePolynomialZp64 = domain.valueOf(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().multiply(el)

    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
    = poly.copy().divide(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: UnivariatePolynomialZp64, el: Long): (UnivariatePolynomialZp64, UnivariatePolynomialZp64)
    = (poly.copy().divide(el), poly.createZero())

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): Long = i.asInstanceOf[Long]
  }

  object UnivariateRingZp64 {
    /**
      * Zp[variable] with specified modulus
      */
    def apply(modulus: Long, variable: String): UnivariateRingZp64 = new UnivariateRingZp64(new IntegersZp64(modulus), variable)

    /**
      * Zp[variable] with specified coefficient ring (Zp)
      */
    def apply(domain: IntegersZp64, variable: String) = new UnivariateRingZp64(domain, variable)
  }

  /**
    * Ring of univariate polynomials over generic domains
    *
    * @param coefficientDomain coefficient ring
    * @param variable          variable
    */
  final case class UnivariateRing[E](coefficientDomain: Ring[E], override val variable: String)
    extends IUnivariateRing[UnivariatePolynomial[E], E](rings.Rings.UnivariateRing(coefficientDomain.domain), variable) {

    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = obj match {
      case poly: UnivariatePolynomial[E] => poly.toString(coefficientDomain, variable)
      case cfx if coefficientDomain.isElement(cfx) => coefficientDomain.show(cfx)
      case _ => super.show(obj)
    }

    /**
      * @inheritdoc
      */
    override def parse(string: String): UnivariatePolynomial[E] = domain.factory().parsePoly(string, coefficientDomain, variable)

    /**
      * @inheritdoc
      */
    override val toString: String = domain.toString(coefficientDomain.toString, variables)

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: E) = domain.factory().createConstant(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().multiply(el)

    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
    = poly.copy().divideExact(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: UnivariatePolynomial[E], el: E): (UnivariatePolynomial[E], UnivariatePolynomial[E])
    = divRem(poly, poly.createConstant(el))

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): E = coefficientDomain.domain.valueOf(i.asInstanceOf[Long])
  }

  private type Ordering = Comparator[DegreeVector[_]]

  /**
    * Ring of multivariate polynomials
    *
    * @param domain    the [[PolynomialRing]]
    * @param variables the variables
    */
  sealed abstract class IMultivariateRing[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (override val domain: poly.PolynomialRing[Poly],
   variables: Array[String],
   ordering: Ordering) extends PolynomialRing[Poly, E](domain, variables)


  /**
    * Zp[variables] with specified modulus
    *
    * @param coefficientDomain coefficient ring
    */
  final case class MultivariateRingZp64
  (coefficientDomain: IntegersZp64, override val variables: Array[String], ordering: Ordering)
    extends IMultivariateRing[MonomialZp64, MultivariatePolynomialZp64, Long](
      rings.Rings.MultivariateRingZp(variables.length, coefficientDomain), variables, ordering) {

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: Long): MultivariatePolynomialZp64 = domain.valueOf(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
    = poly.copy().multiply(el)

    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
    = poly.copy().divide(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: MultivariatePolynomialZp64, el: Long): (MultivariatePolynomialZp64, MultivariatePolynomialZp64)
    = (poly.copy().divide(el), poly.createZero())

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): Long = i.asInstanceOf[Long]
  }

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
    * @param coefficientDomain coefficient ring
    */
  final case class MultivariateRing[E]
  (coefficientDomain: Ring[E], override val variables: Array[String], ordering: Ordering)
    extends IMultivariateRing[Monomial[E], MultivariatePolynomial[E], E](
      rings.Rings.MultivariateRing(variables.length, coefficientDomain), variables, ordering) {
    /**
      * @inheritdoc
      */
    override def show(obj: WithVariables): String = obj match {
      case poly: MultivariatePolynomial[E] => poly.toString(coefficientDomain, variables)
      case el if coefficientDomain.isElement(el) => coefficientDomain.show(el)
      case _ => super.show(obj)
    }

    /**
      * @inheritdoc
      */
    override def parse(string: String): MultivariatePolynomial[E] = domain.factory().parsePoly(string, coefficientDomain, variables)

    /**
      * @inheritdoc
      */
    override val toString: String = domain.toString(coefficientDomain.toString, variables)

    /**
      * Constant polynomial with specified value
      */
    override def getConstant(value: E) = domain.factory().createConstant(value)

    /**
      * Add coefficient ring element
      */
    override def add(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
    = poly.copy().add(el)

    /**
      * Subtract coefficient ring element
      */
    override def subtract(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
    = poly.copy().subtract(el)

    /**
      * Multiply by coefficient ring element
      */
    override def multiply(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
    = poly.copy().multiply(el)

    /**
      * Divide by coefficient ring element
      */
    override def divide(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
    = poly.copy().divideExact(el)

    /**
      * Divide by coefficient ring element
      */
    override def divideAndRemainder(poly: MultivariatePolynomial[E], el: E): (MultivariatePolynomial[E], MultivariatePolynomial[E])
    = divRem(poly, poly.createConstant(el))

    /**
      * Value of integer in coefficient ring
      */
    override def cfValue(i: Int): E = coefficientDomain.valueOf(i.asInstanceOf[Long])
  }

  object MultivariateRing {
    /**
      * Zp[variables] with specified modulus, variables and default ordering (LEX)
      */
    def apply[E](coefficientDomain: Ring[E], variables: Array[String]): MultivariateRing[E]
    = MultivariateRing[E](coefficientDomain, variables, MonomialOrder.LEX)
  }
}
