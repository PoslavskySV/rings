package cc.redberry.rings.scaladsl

import cc.redberry.rings
import cc.redberry.rings._
import cc.redberry.rings.io.{Coder, IParser, IStringifier, Stringifiable}
import cc.redberry.rings.poly.multivar.MonomialOrder

import scala.language.{existentials, implicitConversions, postfixOps}

/**
  * Simple wrapper around [[Ring]] used to unify IPolynomialRing and Ring
  *
  * @param theRing the [[Ring]]
  **/
sealed class Ring[E](val theRing: rings.Ring[E])
  extends Stringifiable[E]
    with IParser[E]
    with RingSupport[E]
    with Serializable {

  override def ringEv(ev: E): Ring[E] = this

  /**
    * Element type
    */
  final type ElementType = E

  /**
    * String from/to conversion for ring elements
    */
  lazy val coder: ACoder[ElementType] = Coder.mkCoder(theRing)

  /**
    * Convert element to string
    */
  final def stringify(e: E): String = coder.stringify(e)

  /**
    * Convert element to string
    */
  final def stringify(e: TraversableOnce[E]): String = e.map(stringify).mkString

  /**
    * Convert element to string
    */
  final def stringify(e: Stringifiable[E]): String = e.toString(coder)

  /**
    * Pretty toString
    */
  @deprecated("use stringify")
  final def show(arg: Any): String =
    if (isElement(arg))
      stringify(arg.asInstanceOf[E])
    else
      arg.toString

  /**
    * @inheritdoc
    **/
  override def toString: String = theRing.toString(coder)

  /**
    * @inheritdoc
    **/
  override def parse(string: String): E = coder parse string

  /** Parse */
  final def apply(string: String): E = parse(string)

  final def apply(int: Int): E = theRing.valueOf(int)

  final def apply(int: BigInt): E = theRing.valueOfBigInteger(int)

  final def apply(int: IntZ): E = theRing.valueOfBigInteger(int)

  final def apply(e: ElementType): E = theRing.valueOf(e)


  final def apply(a: String, b: String): (E, E) = (apply(a), apply(b))

  final def apply(a: String, b: String, c: String): (E, E, E)
  = (apply(a), apply(b), apply(c))

  final def apply(a: String, b: String, c: String, d: String): (E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d))

  final def apply(a: String, b: String, c: String, d: String, e: String): (E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e))

  final def apply(a: String, b: String, c: String, d: String, e: String, f: String): (E, E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e), apply(f))

  final def apply(a: String, b: String, c: String, d: String, e: String, f: String, g: String): (E, E, E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e), apply(f), apply(g))

  final def apply(a: Int, b: Int): (E, E) = (apply(a), apply(b))

  final def apply(a: Int, b: Int, c: Int): (E, E, E)
  = (apply(a), apply(b), apply(c))

  final def apply(a: Int, b: Int, c: Int, d: Int): (E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d))

  final def apply(a: Int, b: Int, c: Int, d: Int, e: Int): (E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e))

  final def apply(a: Int, b: Int, c: Int, d: Int, e: Int, f: Int): (E, E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e), apply(f))

  final def apply(a: Int, b: Int, c: Int, d: Int, e: Int, f: Int, g: Int): (E, E, E, E, E, E, E)
  = (apply(a), apply(b), apply(c), apply(d), apply(e), apply(f), apply(g))

  /** element class  */
  private lazy val eClass: Class[_ <: E] = theRing.getZero.getClass

  /**
    * Reflection: determines whether is element of this
    */
  def isElement(e: Any): Boolean = eClass isAssignableFrom e.getClass

  /**
    * Casts e to element of this
    */
  final def element(e: Any): E = e.asInstanceOf[E]

  final def zero = theRing.getZero
  final def one = theRing.getOne
  final def negativeOne = theRing.getNegativeOne
}

/**
  * Ring of rationals
  */
final case class Frac[E](ring: Ring[E]) extends Ring[Rational[E]](rings.Rings.Frac(scaladsl.ringMethods(ring))) {
  private val rationalsDomain: rings.Rationals[E] = theRing.asInstanceOf[rings.Rationals[E]]

  private[scaladsl] val fracRing: rings.Rationals[E] = theRing.asInstanceOf[rings.Rationals[E]]

  /**
    * String from/to conversion for ring elements
    */
  override lazy val coder: ACoder[ElementType] = Coder.mkRationalsCoder[E](fracRing, ring.coder)

  def apply(a: E): Rational[E] = new Rational[E](rationalsDomain.ring, a)

  def apply(a: E, b: E): (Rational[E], Rational[E]) = (apply(a), apply(b))

  def apply(a: E, b: E, c: E): (Rational[E], Rational[E], Rational[E]) = (apply(a), apply(b), apply(c))

  def apply(a: E, b: E, c: E, d: E): (Rational[E], Rational[E], Rational[E], Rational[E]) = (apply(a), apply(b), apply(c), apply(d))

  def apply(a: E, b: E, c: E, d: E, e: E): (Rational[E], Rational[E], Rational[E], Rational[E], Rational[E]) = (apply(a), apply(b), apply(c), apply(d), apply(e))
}

/**
  * Base class for polynomial rings
  *
  * @param theRing   the [[IPolynomialRing]]
  * @param variables polynomial variables
  * @tparam E coefficient type
  */
abstract class IPolynomialRing[Poly <: IPolynomial[Poly], E]
(override val theRing: rings.poly.IPolynomialRing[Poly],
 val variables: Array[String]) extends Ring[Poly](theRing) {

  /** Set names of variables (new ring will be created) */
  def setVariableNames(newVariables: Array[String]): IPolynomialRing[Poly, E]

  /**
    * Type of coefficients
    */
  final type CoefficientType = E

  /**
    * Type of polynomials
    */
  final type PolyType = ElementType

  /**
    * The coefficient ring
    */
  def cfRing: Ring[E]

  /**
    * Reflection: determines whether is element of this
    */
  override def isElement(e: Any): Boolean
  = super.isElement(e) && e.asInstanceOf[Poly].sameCoefficientRingWith(theRing.getZero)

  /**
    * Index of variable with specified string representation
    */
  final def index(variable: String): Int = variables.indexOf(variable)

  /**
    * Index of variable with specified string representation
    */
  final def variable(variable: String): Int = index(variable)

  /**
    * String representation of i-th variable
    */
  final def variableString(variable: Int): String = variables(variable)

  /**
    * The first variable
    */
  final lazy val `x`: Poly = theRing.variable(0)

  /**
    * The second variable
    */
  final lazy val `y`: Poly = theRing.variable(1)

  /**
    * The third variable
    */
  final lazy val `z`: Poly = theRing.variable(2)

  /**
    * Shortcut for /% operation
    */
  protected[scaladsl] final def divRem(a: Poly, b: Poly): (Poly, Poly) = {
    val qd = theRing.divideAndRemainder(a, b)
    if (qd == null)
      throw new ArithmeticException(s"not divisible with remainder: ${this show a} / ${this show b}")
    (qd(0), qd(1))
  }

  final def apply(value: E): Poly = getConstant(value)

  final def apply(a: E, b: E): (Poly, Poly) = (apply(a), apply(b))

  final def apply(a: E, b: E, c: E): (Poly, Poly, Poly) = (apply(a), apply(b), apply(c))

  final def apply(a: E, b: E, c: E, d: E): (Poly, Poly, Poly, Poly) = (apply(a), apply(b), apply(c), apply(d))

  final def apply(a: E, b: E, c: E, d: E, e: E): (Poly, Poly, Poly, Poly, Poly) = (apply(a), apply(b), apply(c), apply(d), apply(e))

  /**
    * Constant polynomial with specified value
    */
  def getConstant(value: E): Poly

  /**
    * Add coefficient ring element
    */
  def addConstant(poly: Poly, el: E): Poly

  /**
    * Subtract coefficient ring element
    */
  def subtractConstant(poly: Poly, el: E): Poly

  /**
    * Multiply by coefficient ring element
    */
  def multiplyConstant(poly: Poly, el: E): Poly

  /**
    * Divide by coefficient ring element
    */
  def divideConstant(poly: Poly, el: E): Poly

  /**
    * Divide by coefficient ring element
    */
  def divideAndRemainder(poly: Poly, el: E): (Poly, Poly)

  /**
    * Value of integer in coefficient ring
    */
  def cfValue(i: Int): E

  /**
    * Constant coefficient
    */
  def cc(poly: Poly): E

  /**
    * Leading coefficient
    */
  def lc(poly: Poly): E
}

object PolynomialRing {
  def apply[Poly <: IPolynomial[Poly], E](factory: Poly): IPolynomialRing[Poly, E] = factory match {
    case p: UnivariatePolynomialZp64 => UnivariateRingZp64(p.ring, IStringifier.defaultVar()).asInstanceOf[IPolynomialRing[Poly, E]]
    case p: UnivariatePolynomial[E forSome {type E}] => UnivariateRing(p.ring, IStringifier.defaultVar()).asInstanceOf[IPolynomialRing[Poly, E]]
    case p: MultivariatePolynomialZp64 => MultivariateRingZp64(p.ring, IStringifier.defaultVars(p.nVariables), p.ordering).asInstanceOf[IPolynomialRing[Poly, E]]
    case p: MultivariatePolynomial[E forSome {type E}] => MultivariateRing(p.ring, IStringifier.defaultVars(p.nVariables), p.ordering).asInstanceOf[IPolynomialRing[Poly, E]]
    case _ => ???
  }
}

/**
  * Ring of univariate polynomials
  *
  * @param theRing  the [[IPolynomialRing]]
  * @param variable the variable
  */
sealed abstract class IUnivariateRing[Poly <: IUnivariatePolynomial[Poly], E]
(override val theRing: rings.poly.IPolynomialRing[Poly], val variable: String)
  extends IPolynomialRing[Poly, E](theRing, Array(variable)) {

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): IUnivariateRing[Poly, E]

  /**
    * Evaluate poly at a given point
    */
  def eval(poly: Poly, point: E): E

  /**
    * i-th coefficient
    */
  def at(poly: Poly, index: Int): E = cc(poly.getAsPoly(index))

  /**
    * Create univariate polynomial from the array of coefficients
    */
  def create(coefficients: E*): Poly
}

private[scaladsl] sealed abstract class AUnivariateRingZp64
(override val theRing: rings.poly.IPolynomialRing[UnivariatePolynomialZp64], override val variable: String)
  extends IUnivariateRing[UnivariatePolynomialZp64, Long](theRing, variable) {

  /**
    * String from/to conversion for ring elements
    */
  final override lazy val coder: ACoder[ElementType] = Coder.mkUnivariateCoder(theRing, variable)

  /**
    * Constant polynomial with specified value
    */
  override final def getConstant(value: Long): UnivariatePolynomialZp64 = theRing.valueOf(value)

  /**
    * Add coefficient ring element
    */
  override final def addConstant(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
  = theRing valueOf poly.copy().add(el)

  /**
    * Subtract coefficient ring element
    */
  override final def subtractConstant(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
  = theRing valueOf poly.copy().subtract(el)

  /**
    * Multiply by coefficient ring element
    */
  override final def multiplyConstant(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
  = theRing valueOf poly.copy().multiply(el)

  /**
    * Divide by coefficient ring element
    */
  override final def divideConstant(poly: UnivariatePolynomialZp64, el: Long): UnivariatePolynomialZp64
  = theRing valueOf poly.copy().divide(el)

  /**
    * Divide by coefficient ring element
    */
  override final def divideAndRemainder(poly: UnivariatePolynomialZp64, el: Long): (UnivariatePolynomialZp64, UnivariatePolynomialZp64)
  = (divideConstant(poly, el), theRing.getZero)

  /**
    * Value of integer in coefficient ring
    */
  override final def cfValue(i: Int): Long = i.asInstanceOf[Long]

  /**
    * Evaluate poly at a given point
    */
  override final def eval(poly: UnivariatePolynomialZp64, point: Long): Long = poly.evaluate(point)

  /**
    * Constant coefficient
    */
  override final def cc(poly: UnivariatePolynomialZp64): Long = poly.cc()

  /**
    * Leading coefficient
    */
  override final def lc(poly: UnivariatePolynomialZp64): Long = poly.lc()

  /**
    * Create univariate polynomial from the array of coefficients
    */
  override final def create(coefficients: Long*): UnivariatePolynomialZp64 = theRing.valueOf(theRing.factory().createFromArray(coefficients.toArray))
}

private[scaladsl] sealed abstract class AUnivariateRing[E]
(override val theRing: rings.poly.IPolynomialRing[UnivariatePolynomial[E]], override val variable: String)
  extends IUnivariateRing[UnivariatePolynomial[E], E](theRing, variable) {
  /**
    * String from/to conversion for ring elements
    */
  final override lazy val coder: ACoder[ElementType] = Coder.mkUnivariateCoder(theRing, cfRing.coder, variable)

  /**
    * Constant polynomial with specified value
    */
  override final def getConstant(value: E): UnivariatePolynomial[E] = theRing.factory().createConstant(value)

  /**
    * Add coefficient ring element
    */
  override final def addConstant(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
  = theRing valueOf poly.copy().add(el)

  /**
    * Subtract coefficient ring element
    */
  override final def subtractConstant(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
  = theRing valueOf poly.copy().subtract(el)

  /**
    * Multiply by coefficient ring element
    */
  override final def multiplyConstant(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
  = theRing valueOf poly.copy().multiply(el)

  /**
    * Divide by coefficient ring element
    */
  override final def divideConstant(poly: UnivariatePolynomial[E], el: E): UnivariatePolynomial[E]
  = theRing valueOf poly.copy().divideExact(el)

  /**
    * Divide by coefficient ring element
    */
  override final def divideAndRemainder(poly: UnivariatePolynomial[E], el: E): (UnivariatePolynomial[E], UnivariatePolynomial[E])
  = {
    val qd = divRem(poly, poly.createConstant(el))
    (theRing valueOf qd._1, theRing valueOf qd._2)
  }

  /**
    * Value of integer in coefficient ring
    */
  override final def cfValue(i: Int): E = cfRing.valueOf(i.asInstanceOf[Long])

  /**
    * Evaluate poly at a given point
    */
  override final def eval(poly: UnivariatePolynomial[E], point: E): E = poly.evaluate(point)

  /**
    * Constant coefficient
    */
  override final def cc(poly: UnivariatePolynomial[E]): E = poly.cc()

  /**
    * Leading coefficient
    */
  override final def lc(poly: UnivariatePolynomial[E]): E = poly.lc()

  /**
    * Create univariate polynomial from the array of coefficients
    */
  override final def create(coefficients: E*): UnivariatePolynomial[E] = theRing.valueOf(UnivariatePolynomial.apply[E, E](coefficients: _*)(cfRing))
}

sealed trait SimpleFieldExtension[E <: IUnivariatePolynomial[E], C] extends IPolynomialRing[E, C] {
  def implicitConversions: rings.poly.SimpleFieldExtension[E]
}

object SimpleFieldExtension {
  def apply[E <: IUnivariatePolynomial[E], C](field: poly.SimpleFieldExtension[E], variable: String)
  : SimpleFieldExtension[E, C] = JavaConversions.mkScalaFieldExtension(field, variable)

  def apply[E <: IUnivariatePolynomial[E], C](minimalPoly: E, variable: String)
  : SimpleFieldExtension[E, C] = apply(rings.Rings.SimpleFieldExtension(minimalPoly), variable)

  implicit def implicitConversions[E <: IUnivariatePolynomial[E], C](ext: SimpleFieldExtension[E, C]): rings.poly.SimpleFieldExtension[E] =
    ext.implicitConversions
}

/**
  * Galois field with prime base in a range of `(0, 2^63)`
  *
  * @param theRing  the [[rings.poly.FiniteField]]
  * @param variable the variable of univariate polynomials representing this Galois field
  */
final case class GaloisField64
(override val theRing: rings.poly.FiniteField[UnivariatePolynomialZp64], override val variable: String)
  extends AUnivariateRingZp64(theRing, variable) with SimpleFieldExtension[UnivariatePolynomialZp64, Long] {
  /**
    * The coefficient ring
    */
  override val cfRing: Ring[Long] = theRing.factory().ring

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): GaloisField64 = copy(variable = newVariables(0))

  override def implicitConversions: poly.SimpleFieldExtension[UnivariatePolynomialZp64] = theRing
}

object GaloisField64 {
  implicit def asFiniteField(gf: GaloisField64): rings.poly.FiniteField[UnivariatePolynomialZp64] = gf.theRing
}

/**
  * Galois field with arbitrary prime base
  *
  * @param theRing  the [[rings.poly.FiniteField]]
  * @param variable the variable of univariate polynomials representing this Galois field
  */
final case class GaloisField[E]
(override val theRing: rings.poly.FiniteField[UnivariatePolynomial[E]], override val variable: String, _cfRing: Ring[E] = null)
  extends AUnivariateRing[E](theRing, variable) with SimpleFieldExtension[UnivariatePolynomial[E], E] {
  /**
    * The coefficient ring
    */
  override val cfRing: Ring[E] = if (_cfRing == null) theRing.factory().ring else _cfRing

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): GaloisField[E] = copy(variable = newVariables(0))

  override def implicitConversions: poly.SimpleFieldExtension[UnivariatePolynomial[E]] = theRing
}

object GaloisField {
  implicit def asFiniteField[E](gf: GaloisField[E]): rings.poly.FiniteField[UnivariatePolynomial[E]] = gf.theRing
}

object GF {
  /**
    * Create Galois field with cardinality `modulus ^ exponent` represented by univariate polynomials with
    * specified variable
    */
  def apply(modulus: Long, exponent: Int, variable: String = "z")
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
  def apply(modulus: IntZ, exponent: Int, variable: String)
  = GaloisField[IntZ](rings.Rings.GF(modulus, exponent), variable)

  /**
    * Create Galois field from the specified irreducible polynomial and represented by univariate polynomials with
    * specified variable
    */
  def apply[E](irreducible: UnivariatePolynomial[E], variable: String)(implicit cfRing: Ring[E])
  = GaloisField[E](rings.Rings.GF(irreducible), variable, _cfRing = cfRing)
}

/**
  * Ring of Zp[x] polynomials
  *
  * @param cfRingZp64 coefficient ring
  * @param variable   variable
  */
final case class UnivariateRingZp64(cfRingZp64: IntegersZp64, override val variable: String)
  extends AUnivariateRingZp64(rings.Rings.UnivariateRingZp64(cfRingZp64), variable) {
  val modulus: Long = cfRingZp64.modulus
  /**
    * The coefficient ring
    */
  override val cfRing: Ring[Long] = cfRingZp64

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): UnivariateRingZp64 = copy(variable = newVariables(0))
}

object UnivariateRingZp64 {
  /**
    * Zp[variable] with specified modulus
    */
  def apply(modulus: Long, variable: String): UnivariateRingZp64 = new UnivariateRingZp64(new IntegersZp64(modulus), variable)
}

/**
  * Ring of univariate polynomials over generic domains
  *
  * @param cfRing   coefficient ring
  * @param variable variable
  */
final case class UnivariateRing[E](override val cfRing: Ring[E], override val variable: String)
  extends AUnivariateRing[E](rings.Rings.UnivariateRing(cfRing.theRing), variable) {
  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): UnivariateRing[E] = copy(variable = newVariables(0))
}

/**
  * Algebraic number field represented as simple field extension
  *
  * @param theRing  the [[rings.poly.FiniteField]]
  * @param variable the variable that represent extension generator
  */
final case class AlgebraicNumberField[E]
(override val theRing: rings.poly.AlgebraicNumberField[UnivariatePolynomial[E]], override val variable: String, _cfRing: Ring[E] = null)
  extends AUnivariateRing[E](theRing, variable) with SimpleFieldExtension[UnivariatePolynomial[E], E] {
  /**
    * The coefficient ring
    */
  override val cfRing: Ring[E] = if (_cfRing == null) theRing.factory().ring else _cfRing

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): AlgebraicNumberField[E] = copy(variable = newVariables(0))

  override def implicitConversions: poly.SimpleFieldExtension[UnivariatePolynomial[E]] = theRing
}

object AlgebraicNumberField {
  def apply[E](minimalPolynomial: UnivariatePolynomial[E], variable: String): AlgebraicNumberField[E]
  = AlgebraicNumberField(rings.Rings.AlgebraicNumberField(minimalPolynomial), variable)

  implicit def asAlgebraicNumberField[E](field: AlgebraicNumberField[E]): rings.poly.AlgebraicNumberField[UnivariatePolynomial[E]] = field.theRing
}

/**
  * Ring of multivariate polynomials
  *
  * @param theRing   the [[IPolynomialRing]]
  * @param variables the variables
  */
sealed abstract class IMultivariateRing[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
(override val theRing: rings.poly.IPolynomialRing[Poly], variables: Array[String],
 val ordering: Ordering) extends IPolynomialRing[Poly, E](theRing, variables) {

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): IMultivariateRing[Term, Poly, E]

  /**
    * The type of monomials
    */
  type MonomialType = Term

  /**
    * Evaluate poly for given variable
    */
  def eval(poly: Poly, variable: Int, value: E): Poly

  //    /**
  //      * Evaluate poly for given variables
  //      */
  //    def eval(poly: Poly, variables: Array[Int], values: Array[E]): Poly
}


/**
  * Zp[variables] with specified modulus
  *
  * @param coefficientDomain coefficient ring
  */
final case class MultivariateRingZp64
(coefficientDomain: IntegersZp64, override val variables: Array[String], override val ordering: Ordering)
  extends IMultivariateRing[MonomialZp64, MultivariatePolynomialZp64, Long](
    rings.Rings.MultivariateRingZp64(variables.length, coefficientDomain, ordering), variables, ordering) {

  private[scaladsl] val multivariateRing: poly.MultivariateRing[PolyType] = theRing.asInstanceOf[rings.poly.MultivariateRing[PolyType]]

  /**
    * String from/to conversion for ring elements
    */
  override lazy val coder: ACoder[ElementType] = Coder.mkMultivariateCoder[MonomialType, PolyType](multivariateRing, variables: _*)

  /**
    * The coefficient ring
    */
  override val cfRing: Ring[Long] = coefficientDomain

  /**
    * Constant polynomial with specified value
    */
  override def getConstant(value: Long): MultivariatePolynomialZp64 = theRing.valueOf(value)

  /**
    * Add coefficient ring element
    */
  override def addConstant(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
  = poly.copy().add(el)

  /**
    * Subtract coefficient ring element
    */
  override def subtractConstant(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
  = poly.copy().subtract(el)

  /**
    * Multiply by coefficient ring element
    */
  override def multiplyConstant(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
  = poly.copy().multiply(el)

  /**
    * Divide by coefficient ring element
    */
  override def divideConstant(poly: MultivariatePolynomialZp64, el: Long): MultivariatePolynomialZp64
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

  /**
    * Evaluate poly for given variable
    */
  override def eval(poly: MultivariatePolynomialZp64, variable: Int, value: Long): MultivariatePolynomialZp64
  = poly.evaluate(variable, value)

  /**
    * Constant coefficient
    */
  override def cc(poly: MultivariatePolynomialZp64): Long = poly.cc()

  /**
    * Leading coefficient
    */
  override def lc(poly: MultivariatePolynomialZp64): Long = poly.lc()

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): MultivariateRingZp64 = copy(variables = newVariables)
}

object MultivariateRingZp64 {
  /**
    * Zp[variables] with specified modulus, variables and ordering
    */
  def apply(modulus: Long, variables: Array[String], ordering: Ordering): MultivariateRingZp64
  = MultivariateRingZp64(new IntegersZp64(modulus), variables, ordering)

  /**
    * Zp[variables] with specified modulus, variables and default ordering (MonomialOrder.DEFAULT)
    */
  def apply(modulus: Long, variables: Array[String]): MultivariateRingZp64
  = MultivariateRingZp64(new IntegersZp64(modulus), variables, MonomialOrder.DEFAULT)
}

/**
  * Ring of multivariate polynomials over generic domains
  *
  * @param coefficientDomain coefficient ring
  */
final case class MultivariateRing[E]
(coefficientDomain: Ring[E], override val variables: Array[String], override val ordering: Ordering)
  extends IMultivariateRing[Monomial[E], MultivariatePolynomial[E], E](
    rings.Rings.MultivariateRing(variables.length, scaladsl.ringMethods(coefficientDomain), ordering), variables, ordering) {
  private[scaladsl] val multivariateRing: poly.MultivariateRing[PolyType] = theRing.asInstanceOf[rings.poly.MultivariateRing[PolyType]]

  /**
    * String from/to conversion for ring elements
    */
  override lazy val coder: ACoder[ElementType] = Coder.mkMultivariateCoder[E](multivariateRing, cfRing.coder, variables: _*)

  /**
    * The coefficient ring
    */
  override val cfRing: Ring[E] = coefficientDomain

  /**
    * Constant polynomial with specified value
    */
  override def getConstant(value: E) = theRing.factory().createConstant(value)

  /**
    * Add coefficient ring element
    */
  override def addConstant(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
  = poly.copy().add(el)

  /**
    * Subtract coefficient ring element
    */
  override def subtractConstant(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
  = poly.copy().subtract(el)

  /**
    * Multiply by coefficient ring element
    */
  override def multiplyConstant(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
  = poly.copy().multiply(el)

  /**
    * Divide by coefficient ring element
    */
  override def divideConstant(poly: MultivariatePolynomial[E], el: E): MultivariatePolynomial[E]
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

  /**
    * Evaluate poly for given variable
    */
  override def eval(poly: MultivariatePolynomial[E], variable: Int, value: E): MultivariatePolynomial[E]
  = poly.evaluate(variable, value)

  /**
    * Constant coefficient
    */
  override def cc(poly: MultivariatePolynomial[E]): E = poly.cc()

  /**
    * Leading coefficient
    */
  override def lc(poly: MultivariatePolynomial[E]): E = poly.cc()

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): MultivariateRing[E] = copy(variables = newVariables)
}

object MultivariateRing {
  /**
    * Zp[variables] with specified modulus, variables and default ordering (MonomialOrder.DEFAULT)
    */
  def apply[E](coefficientDomain: Ring[E], variables: Array[String]): MultivariateRing[E]
  = MultivariateRing[E](coefficientDomain, variables, MonomialOrder.DEFAULT)

  /**
    * Zp[variables] with specified modulus, variables and default ordering (MonomialOrder.DEFAULT)
    */
  def apply[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](factory: Poly): IMultivariateRing[Term, Poly, E]
  = factory match {
    case _: poly.multivar.MultivariatePolynomialZp64 => {
      val f = factory.asInstanceOf[MultivariatePolynomialZp64]
      new MultivariateRingZp64(f.ring, IStringifier.defaultVars(f.nVariables), f.ordering).asInstanceOf[IMultivariateRing[Term, Poly, E]]
    }
    case _: poly.multivar.MultivariatePolynomial[_] => {
      val f = factory.asInstanceOf[MultivariatePolynomial[_]]
      new MultivariateRing[E](asRing(f.ring).asInstanceOf[Ring[E]], IStringifier.defaultVars(f.nVariables), f.ordering).asInstanceOf[IMultivariateRing[Term, Poly, E]]
    }
  }
}

/**
  * Ideal in multivariate polynomial ring
  */
final case class Ideal[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
(ring: IMultivariateRing[Term, Poly, E], theIdeal: poly.multivar.Ideal[Term, Poly]) {
  /**
    * Set the monomial order used for Groebner basis of this ideal
    */
  def changeOrder(newOrder: Ordering) = Ideal(ring, theIdeal.changeOrder(newOrder))

  def +(oth: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] = union(oth)

  def +(oth: Poly): Ideal[Term, Poly, E] = union(oth)

  def *(oth: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] = multiply(oth)

  def *(oth: Poly): Ideal[Term, Poly, E] = multiply(oth)

  def :/(oth: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] = quotient(oth)

  def :/(oth: Poly): Ideal[Term, Poly, E] = quotient(oth)

  def ∪(oth: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] = union(oth)

  def ∪(oth: Poly): Ideal[Term, Poly, E] = union(oth)

  def ∩(oth: Ideal[Term, Poly, E]): Ideal[Term, Poly, E] = intersection(oth)

  lazy val groebnerBasis: Seq[Poly] = {
    import scala.collection.JavaConverters._
    theIdeal.getGroebnerBasis.asScala.toList
  }

  /**
    * Returns the union of this and oth
    */
  def union(oth: Ideal[Term, Poly, E]) = Ideal(ring, theIdeal.union(oth))

  /**
    * Returns the union of this and oth
    */
  def union(oth: Poly) = Ideal(ring, theIdeal.union(oth))

  /**
    * Returns the intersection of this and oth
    */
  def intersection(oth: Ideal[Term, Poly, E]) = Ideal(ring, theIdeal.intersection(oth))

  /**
    * Returns the product of this and oth
    */
  def multiply(oth: Ideal[Term, Poly, E]) = Ideal(ring, theIdeal.multiply(oth))

  /**
    * Returns the product of this and oth
    */
  def multiply(oth: Poly) = Ideal(ring, theIdeal.multiply(oth))

  /**
    * Returns this in a power of exponent
    */
  def pow(exponent: Integer) = Ideal(ring, theIdeal.pow(exponent))

  /**
    * Returns the quotient this : oth
    */
  def quotient(oth: Ideal[Term, Poly, E]) = Ideal(ring, theIdeal.quotient(oth))

  /**
    * Returns the quotient this : oth
    */
  def quotient(oth: Poly) = Ideal(ring, theIdeal.quotient(oth))

  /**
    * Ideal of leading terms
    */
  lazy val ltIdeal: Ideal[Term, Poly, E] = Ideal(ring, theIdeal.ltIdeal())

  override def toString = theIdeal.toString(ring.coder)
}

object Ideal {

  import scala.collection.JavaConverters._

  def apply[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (ring: IMultivariateRing[Term, Poly, E], generators: Seq[Poly], monomialOrder: Ordering)
  : Ideal[Term, Poly, E] =
    new Ideal[Term, Poly, E](ring, poly.multivar.Ideal.create[Term, Poly](generators.toList.asJava, monomialOrder))

  def apply[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (ring: IMultivariateRing[Term, Poly, E], generators: Seq[Poly])
  : Ideal[Term, Poly, E] = apply(ring, generators, ring.factory().ordering)

  def apply[E](generators: Seq[MultivariatePolynomial[E]], monomialOrder: Ordering = null)
              (implicit ring: IMultivariateRing[Monomial[E], MultivariatePolynomial[E], E] = null)
  : Ideal[Monomial[E], MultivariatePolynomial[E], E] = {
    val r: IMultivariateRing[Monomial[E], MultivariatePolynomial[E], E] =
      if (ring != null)
        ring
      else
        MultivariateRing(generators(0))

    val m =
      if (monomialOrder == null)
        r.theRing.factory().ordering
      else monomialOrder

    new Ideal[Monomial[E], MultivariatePolynomial[E], E](r,
      poly.multivar.Ideal.create[Monomial[E], MultivariatePolynomial[E]](generators.toList.asJava, m))
  }
}

object IdealZp64 {

  import scala.collection.JavaConverters._

  def apply(generators: Seq[MultivariatePolynomialZp64], monomialOrder: Ordering = null)
           (implicit ring: IMultivariateRing[MonomialZp64, MultivariatePolynomialZp64, Long] = null)
  : Ideal[MonomialZp64, MultivariatePolynomialZp64, Long] = {
    val r: IMultivariateRing[MonomialZp64, MultivariatePolynomialZp64, Long] =
      if (ring != null)
        ring
      else
        MultivariateRing(generators(0))

    val m =
      if (monomialOrder == null)
        r.theRing.factory().ordering
      else monomialOrder

    new Ideal[MonomialZp64, MultivariatePolynomialZp64, Long](r,
      poly.multivar.Ideal.create[MonomialZp64, MultivariatePolynomialZp64](generators.toList.asJava, m))
  }
}

private[scaladsl]
abstract class MultivariateRingWrapper[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
(baseRing: rings.poly.IPolynomialRing[Poly], variables: Array[String])
  extends IMultivariateRing[Term, Poly, E](baseRing, variables, baseRing.factory().ordering) {
  protected[scaladsl] lazy val helperRing: IMultivariateRing[Term, Poly, E] = MultivariateRing.apply(baseRing.factory())

  /**
    * Evaluate poly for given variable
    */
  override def eval(poly: Poly, variable: Int, value: E) = helperRing.eval(poly, variable, value)

  /**
    * The coefficient ring
    */
  override def cfRing = helperRing.cfRing

  /**
    * Constant polynomial with specified value
    */
  override def getConstant(value: E) = helperRing.getConstant(value)

  /**
    * Add coefficient ring element
    */
  override def addConstant(poly: Poly, el: E) = helperRing.addConstant(poly, el)

  /**
    * Subtract coefficient ring element
    */
  override def subtractConstant(poly: Poly, el: E) = helperRing.subtractConstant(poly, el)

  /**
    * Multiply by coefficient ring element
    */
  override def multiplyConstant(poly: Poly, el: E) = helperRing.multiplyConstant(poly, el)

  /**
    * Divide by coefficient ring element
    */
  override def divideConstant(poly: Poly, el: E) = helperRing.divideConstant(poly, el)

  /**
    * Divide by coefficient ring element
    */
  override def divideAndRemainder(poly: Poly, el: E) = helperRing.divideAndRemainder(poly, el)

  /**
    * Value of integer in coefficient ring
    */
  override def cfValue(i: Int) = helperRing.cfValue(i)

  /**
    * Constant coefficient
    */
  override def cc(poly: Poly) = helperRing.cc(poly)

  /**
    * Leading coefficient
    */
  override def lc(poly: Poly) = helperRing.lc(poly)
}

/**
  * Multivariate quotient ring
  */
final case class QuotientRing[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
(baseRing: IMultivariateRing[Term, Poly, E], ideal: Ideal[Term, Poly, E])
  extends MultivariateRingWrapper[Term, Poly, E](
    rings.Rings.QuotientRing[Term, Poly](baseRing.theRing.asInstanceOf[rings.poly.MultivariateRing[Poly]], ideal.theIdeal),
    baseRing.variables) {
  /**
    * String from/to conversion for ring elements
    */
  override lazy val coder: ACoder[ElementType] = baseRing.coder

  /**
    * @inheritdoc
    **/
  override def parse(string: String): Poly = apply(coder.parse(string))

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String]): QuotientRing[Term, Poly, E]
  = copy(baseRing.setVariableNames(newVariables))
}

final case class MultipleFieldExtension[
Term <: AMonomial[Term],
mPoly <: AMultivariatePolynomial[Term, mPoly],
sPoly <: IUnivariatePolynomial[sPoly],
E](override val theRing: rings.poly.MultipleFieldExtension[Term, mPoly, sPoly],
   override val variables: Array[String])
  extends MultivariateRingWrapper[Term, mPoly, E](theRing, variables) {

  /**
    * String from/to conversion for ring elements
    */
  override lazy val coder: ACoder[mPoly] = Coder.mkMultipleExtensionCoder[Term, mPoly, sPoly](theRing, variables: _*)

  /**
    * @inheritdoc
    **/
  override def parse(string: String): mPoly = coder.parse(string)

  /**
    * Returns the isomorphic simple field extension
    */
  def getSimpleExtension(variable: String = "gamma")
  : SimpleFieldExtension[sPoly, E] = {
    val r: SimpleFieldExtension[sPoly, E] = SimpleFieldExtension(theRing.getSimpleExtension, variable)
    for (i <- 0 until theRing.nVariables()) {
      r.coder.bindAlias(variableString(i), theRing.getGeneratorRep(i))
    }
    r
  }

  /**
    * Adds algebraic element given by its minimal polynomial (not checked that it is irreducible) to this.
    */
  def joinAlgebraicElement(algebraicElement: UnivariatePolynomial[mPoly], variable: String)
  : MultipleFieldExtension[Term, mPoly, sPoly, E] = {
    val r: MultipleFieldExtension[Term, mPoly, sPoly, E] = MultipleFieldExtension(theRing.joinAlgebraicElement(algebraicElement), variables :+ variable)
    r.coder.withEncoder(this.coder)
    //r.coder.getBindings.putAll(this.coder.getBindings)
    r
  }

  /**
    * Adds algebraic element given by its minimal polynomial (not checked that it is irreducible) to this.
    */
  def joinAlgebraicElement(algebraicElement: sPoly, variable: String)
  : MultipleFieldExtension[Term, mPoly, sPoly, E] = {
    val r: MultipleFieldExtension[Term, mPoly, sPoly, E] = MultipleFieldExtension(theRing.joinAlgebraicElement(algebraicElement), variables :+ variable)
    r.coder.withEncoder(this.coder)
    //r.coder.getBindings.putAll(this.coder.getBindings)
    r
  }

  /** Set names of variables (new ring will be created) */
  override def setVariableNames(newVariables: Array[String])
  : MultipleFieldExtension[Term, mPoly, sPoly, E] = copy(variables = newVariables)
}

object MultipleFieldExtension {
  def fromSimpleExtension[Term <: AMonomial[Term],
  mPoly <: AMultivariatePolynomial[Term, mPoly],
  uPoly <: IUnivariatePolynomial[uPoly],
  E](simpleExt: SimpleFieldExtension[uPoly, E])
  : MultipleFieldExtension[Term, mPoly, uPoly, E] = JavaConversions.mkSimpleToMultiple(simpleExt)

  def apply(simpleExt: GaloisField64)
  : MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64, Long] = {
    val mExt: poly.MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64] =
      simpleExt
        .implicitConversions
        .asMultipleExtension()

    MultipleFieldExtension(mExt, simpleExt.variables)
  }

  def apply[E](simpleExt: GaloisField[E])
  : MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E], E] = {
    val mExt: poly.MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E]] =
      simpleExt
        .implicitConversions
        .asMultipleExtension()

    MultipleFieldExtension(mExt, simpleExt.variables)
  }

  def apply[E](simpleExt: AlgebraicNumberField[E])
  : MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E], E] = {
    val mExt: poly.MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E]] =
      simpleExt
        .implicitConversions
        .asMultipleExtension()

    MultipleFieldExtension(mExt, simpleExt.variables)
  }

  def apply[E](generators: Array[UnivariatePolynomial[E]],
               variables: Array[String]):
  MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E], E] =
    MultipleFieldExtension(rings.Rings.MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E]](generators: _*), variables)

  def apply(generators: Array[UnivariatePolynomialZp64],
            variables: Array[String]): MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64, Long] =
    MultipleFieldExtension(rings.Rings.MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64](generators: _*), variables)

  def apply[E](generator: UnivariatePolynomial[E],
               variable: String):
  MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E], E] =
    MultipleFieldExtension(rings.Rings.MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E]](generator), Array(variable))

  def apply(generator: UnivariatePolynomialZp64,
            variable: String): MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64, Long] =
    MultipleFieldExtension(rings.Rings.MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64](generator), Array(variable))


  implicit def implicitConversions[
  Term <: AMonomial[Term],
  mPoly <: AMultivariatePolynomial[Term, mPoly],
  uPoly <: IUnivariatePolynomial[uPoly],
  E](multipleFieldExtension: MultipleFieldExtension[Term, mPoly, uPoly, E])
  : rings.poly.MultipleFieldExtension[Term, mPoly, uPoly] = multipleFieldExtension.theRing
}
