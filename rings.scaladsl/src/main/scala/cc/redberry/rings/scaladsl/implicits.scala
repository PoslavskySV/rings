package cc.redberry.rings.scaladsl

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.IPolynomial
import cc.redberry.rings.poly.multivar._
import cc.redberry.rings.poly.univar.{IUnivariatePolynomial, UnivariatePolynomial, UnivariatePolynomialZp64}
import cc.redberry.rings.{IntegersZp64, Rational}

import scala.language.{existentials, implicitConversions}

/**
  *
  *
  * @since 1.0
  */
object implicits {
  type Integer = BigInteger

  import Rings._

  sealed class RingOperators[E](self: E)(implicit ring: Ring[E]) {
    protected[implicits] def constant(value: Int): E = ring.valueOf(value)

    protected[implicits] def constant(value: Long): E = ring.valueOf(value)

    def +(other: Int): E = ring.add(self, ring.valueOf(other))

    def +(other: E): E = ring.add(self, other)

    def +=(other: E): E = ring.addMutable(self, other)

    def -(other: E): E = ring.subtract(self, other)

    def -=(other: E): E = ring.subtractMutable(self, other)

    def *(other: E): E = ring.multiply(self, other)

    def *=(other: E): E = ring.multiplyMutable(self, other)

    def **(exponent: Int): E = ring.pow(self, exponent)

    def ^(exponent: Int): E = this.**(exponent)

    def unary_- : E = ring.negate(self)

    def unary_-= : E = ring.negateMutable(self)

    def unary_+ : E = ring.copy(self)

    def ++ : E = ring.increment(self)

    def -- : E = ring.decrement(self)

    def gcd(other: E): E = ring.gcd(self, other)

    def <(other: E): Boolean = ring.compare(self, other) < 0

    def <=(other: E): Boolean = ring.compare(self, other) <= 0

    def >(other: E): Boolean = ring.compare(self, other) > 0

    def >=(other: E): Boolean = ring.compare(self, other) >= 0

    def /%(other: E): (E, E) = {
      val qd = ring.divideAndRemainder(self, other)
      (qd(0), qd(1))
    }

    def %(other: E): E = ring.remainder(self, other)

    def /(other: E): E = ring.divideExact(self, other)

    def ===(other: Int): Boolean = self == constant(other)

    def ===(other: Long): Boolean = self == constant(other)

    def ===(other: E): Boolean = self == other
  }

  implicit def withRationalOperators[E]
  (self: Rational[E])(implicit ring: Rationals[E]): RingOperators[Rational[E]] = new RingOperators(self)(ring)

  @specialized(Long)
  private final case class
  RingOperatorsZp64(self: Long, domain: IntegersZp64) extends RingOperators[Long](domain.modulus(self))(null) {
    override protected[implicits] def constant(value: Int): Long = constant(value.asInstanceOf[Long])

    override protected[implicits] def constant(value: Long): Long = domain.modulus(value)

    override def +(other: Long): Long = domain.add(self, other)

    override def -(other: Long): Long = domain.subtract(self, other)

    override def *(other: Long): Long = domain.multiply(self, other)

    override def **(exponent: Int): Long = domain.powMod(self, exponent)

    override def ^(exponent: Int): Long = domain.powMod(self, exponent)

    override def unary_- : Long = domain.negate(self)

    override def unary_+ : Long = self

    override def ++ : Long = domain.add(self, 1)

    override def -- : Long = domain.subtract(self, 1)

    override def <(other: Long): Boolean = self < domain.modulus(other)

    override def <=(other: Long): Boolean = self <= domain.modulus(other)

    override def >(other: Long): Boolean = self > domain.modulus(other)

    override def >=(other: Long): Boolean = self >= domain.modulus(other)

    override def /%(other: Long): (Long, Long) = (domain.divide(self, other), 0)

    override def /(other: Long): Long = domain.divide(self, other)

    override def ===(other: Int): Boolean = self == domain.modulus(other)

    override def ===(other: Long): Boolean = self == domain.modulus(other)
  }

  sealed abstract class IPolynomialOperators[Poly <: IPolynomial[Poly], E]
  (self: Poly)(implicit ring: PolynomialRing[Poly, E]) extends RingOperators[Poly](self)(ring) {

    def show(): String = ring.show(self)

    def +(element: E): Poly = ring.add(self, element)

    def -(element: E): Poly = ring.subtract(self, element)

    def *(element: E): Poly = ring.multiply(self, element)

    def /(element: E): Poly = ring.divide(self, element)

    def :=(other: Poly): Unit = self.set(other)

    def /%(element: E)(implicit dummy: DummyImplicit): (Poly, Poly) = ring.divideAndRemainder(self, element)

    def lc(): E

    def cc(): E
  }

  sealed abstract class IUnivariateOperators[Poly <: IUnivariatePolynomial[Poly], E](self: Poly)(implicit ring: PolynomialRing[Poly, E])
    extends IPolynomialOperators[Poly, E](self)(ring) {

    def <<(offset: Int): Poly = ring.valueOf(self.copy().shiftLeft(offset))

    def >>(offset: Int): Poly = ring.valueOf(self.copy().shiftRight(offset))

    def @@(index: Int): E = apply(index)

    def apply(from: Int, to: Int): Poly = ring.valueOf(self.copy().getRange(from, to))

    def apply(poly: Poly): Poly = ring.valueOf(self.composition(poly))

    def apply(index: Int): E

    def eval(point: E): E
  }

  final case class UnivariateZp64Operators
  (self: UnivariatePolynomialZp64)(implicit ring: IUnivariateRing[UnivariatePolynomialZp64, Long])
    extends IUnivariateOperators[UnivariatePolynomialZp64, Long](self)(ring) {

    override def lc(): Long = self.lc()

    override def cc(): Long = self.cc()

    override def apply(index: Int): Long = self.get(index)

    override def eval(point: Long): Long = self.evaluate(point)
  }

  final case class UnivariateOperators[E]
  (self: UnivariatePolynomial[E])(implicit ring: IUnivariateRing[UnivariatePolynomial[E], E])
    extends IUnivariateOperators[UnivariatePolynomial[E], E](self)(ring) {

    override def lc(): E = self.lc()

    override def cc(): E = self.cc()

    override def apply(index: Int): E = self.get(index)

    override def eval(point: E): E = self.evaluate(point)
  }

  private[implicits] sealed abstract class IMultivariateOperators[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
  (self: Poly)(implicit ring: PolynomialRing[Poly, E])
    extends IPolynomialOperators[Poly, E](self)(ring) {
    def +(other: Term): Poly = ring.theRing.valueOf(self.copy().add(other))

    def -(other: Term): Poly = ring.theRing.valueOf(self.copy().subtract(other))

    def *(other: Term): Poly = ring.theRing.valueOf(self.copy().multiply(other))

    def /(other: Term): Poly = {
      val r = self.copy().divideOrNull(other)
      if (r == null)
        throw new ArithmeticException(s"not divisible: $self / $other")
      else
        ring.theRing.valueOf(r)
    }

    def swapVariables(i: Int, j: Int): Poly

    def swapVariables(tup: (String, String)): Poly = swapVariables(ring.index(tup._1), ring.index(tup._2))

    def eval(variable: String, value: E): Poly

    final def eval(substitution: (String, E))(implicit dummyImplicit: DummyImplicit): Poly
    = eval(substitution._1, substitution._2)

    final def eval(substitution: (String, Int)): Poly
    = eval((substitution._1, ring.cfValue(substitution._2)))

    final def apply(substitution: (String, E))(implicit dummyImplicit: DummyImplicit): Poly = eval(substitution)

    final def apply(substitution: (String, Int)): Poly
    = eval(substitution)
  }

  final case class MultivariateZp64Operators
  (self: MultivariatePolynomialZp64)(implicit ring: PolynomialRing[MultivariatePolynomialZp64, Long])
    extends IMultivariateOperators[MonomialZp64, MultivariatePolynomialZp64, Long](self)(ring) {

    override def lc(): Long = self.lc()

    override def cc(): Long = self.cc()

    override def swapVariables(i: Int, j: Int): MultivariatePolynomialZp64 = AMultivariatePolynomial.swapVariables[MonomialZp64, MultivariatePolynomialZp64](self, i, j)

    override def eval(variable: String, value: Long): MultivariatePolynomialZp64 = self.evaluate(ring.index(variable), value)

    def eval(substitutions: (String, Long)*): MultivariatePolynomialZp64 = self.evaluate(substitutions.map(t => ring.index(t._1)).toArray, substitutions.map(_._2).toArray)
  }

  final case class MultivariateOperators[E]
  (self: MultivariatePolynomial[E])(implicit ring: PolynomialRing[MultivariatePolynomial[E], E])
    extends IMultivariateOperators[Monomial[E], MultivariatePolynomial[E], E](self)(ring) {

    override def lc(): E = self.lc()

    override def cc(): E = self.cc()

    override def swapVariables(i: Int, j: Int): MultivariatePolynomial[E] = AMultivariatePolynomial.swapVariables[Monomial[E], MultivariatePolynomial[E]](self, i, j)

    override def eval(variable: String, value: E): MultivariatePolynomial[E] = self.evaluate(ring.index(variable), value)

  }

  final case class StringSwap(self: String) {
    def <>(oth: String): (String, String) = (self, oth)

    def ->(oth: Long): (String, Long) = (self, oth)

    def ->(oth: Int): (String, Int) = (self, oth)

    def ->[E](oth: E): (String, E) = (self, oth)
  }

  implicit def withStringSwap(self: String): StringSwap = StringSwap(self)

  object univariate {
    implicit def withUnivariateOperators[Poly <: IUnivariatePolynomial[Poly], E]
    (self: Poly)(implicit ring: IUnivariateRing[Poly, E]): IUnivariateOperators[Poly, E]
    = self match {
      case p: UnivariatePolynomialZp64 => UnivariateZp64Operators(p)(ring.asInstanceOf[IUnivariateRing[UnivariatePolynomialZp64, Long]]).asInstanceOf[IUnivariateOperators[Poly, E]]
      case p: UnivariatePolynomial[E] => UnivariateOperators[E](p)(ring.asInstanceOf[IUnivariateRing[UnivariatePolynomial[E], E]]).asInstanceOf[IUnivariateOperators[Poly, E]]
      case _ => throw new RuntimeException("unknown type: " + self.getClass + "  " + self.coefficientRingToString())
    }
  }

  object multivariate {
    implicit def withMultivariateOperators[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly], E]
    (self: Poly)(implicit ring: IMultivariateRing[Term, Poly, E]): IMultivariateOperators[Term, Poly, E]
    = self match {
      case p: MultivariatePolynomialZp64 => MultivariateZp64Operators(p)(ring.asInstanceOf[PolynomialRing[MultivariatePolynomialZp64, Long]]).asInstanceOf[IMultivariateOperators[Term, Poly, E]]
      case p: MultivariatePolynomial[E] => MultivariateOperators[E](p)(ring.asInstanceOf[PolynomialRing[MultivariatePolynomial[E], E]]).asInstanceOf[IMultivariateOperators[Term, Poly, E]]
      case _ => throw new RuntimeException("unknown type: " + self.getClass + "  " + self.coefficientRingToString())
    }
  }

  object generic {
    object rings {
      implicit def withRingOperators[E]
      (self: E)(implicit ring: Ring[E]): RingOperators[E] = new RingOperators(self)(ring)
    }
    object polynomials {
      implicit def withPolynomialOperators[Poly <: IPolynomial[Poly], E]
      (self: Poly)(implicit ring: PolynomialRing[Poly, E]): IPolynomialOperators[Poly, E]
      = self match {
        case u64: UnivariatePolynomialZp64 => UnivariateZp64Operators(u64)(ring.asInstanceOf[IUnivariateRing[UnivariatePolynomialZp64, Long]]).asInstanceOf[IPolynomialOperators[Poly, E]]
        case ue: UnivariatePolynomial[E] => UnivariateOperators[E](ue)(ring.asInstanceOf[IUnivariateRing[UnivariatePolynomial[E], E]]).asInstanceOf[IPolynomialOperators[Poly, E]]
        case m64: MultivariatePolynomialZp64 => MultivariateZp64Operators(m64)(ring.asInstanceOf[PolynomialRing[MultivariatePolynomialZp64, Long]]).asInstanceOf[IPolynomialOperators[Poly, E]]
        case me: MultivariatePolynomial[E] => MultivariateOperators[E](me)(ring.asInstanceOf[PolynomialRing[MultivariatePolynomial[E], E]]).asInstanceOf[IPolynomialOperators[Poly, E]]
        case _ => throw new RuntimeException("unknown type: " + self.getClass + "  " + self.coefficientRingToString())
      }
    }
  }

  implicit def withUnivariateOperatorsZp64
  (poly: UnivariatePolynomialZp64)
  (implicit ring: UnivariateRingZp64) =
    UnivariateZp64Operators(poly)(ring)

  implicit def withUnivariateOperatorsE[E]
  (poly: UnivariatePolynomial[E])
  (implicit ring: UnivariateRing[E]) =
    UnivariateOperators[E](poly)(ring)

  implicit def withMultivariateOperatorsZp64
  (poly: MultivariatePolynomialZp64)
  (implicit ring: MultivariateRingZp64) =
    MultivariateZp64Operators(poly)(ring)

  implicit def withMultivariateOperatorsE[E]
  (poly: MultivariatePolynomial[E])
  (implicit ring: MultivariateRing[E]) =
    MultivariateOperators(poly)(ring)

  implicit def withRingOperators[Poly <: IPolynomial[Poly], E]
  (self: E)(implicit ring: PolynomialRing[Poly, E]): RingOperators[E] =
    ring match {
      case r: GaloisField64 => RingOperatorsZp64(self.asInstanceOf[Long], r.theRing.factory.ring).asInstanceOf[RingOperators[E]]
      case r: GaloisField[E] => new RingOperators[E](self)(ringAsScalaRing(r.factory.ring))
      case r: UnivariateRingZp64 => RingOperatorsZp64(self.asInstanceOf[Long], r.coefficientDomain).asInstanceOf[RingOperators[E]]
      case r: UnivariateRing[E] => new RingOperators(self)(r.coefficientDomain)
      case r: MultivariateRingZp64 => RingOperatorsZp64(self.asInstanceOf[Long], r.coefficientDomain).asInstanceOf[RingOperators[E]]
      case r: MultivariateRing[E] => new RingOperators(self)(r.coefficientDomain)
    }

  final case class NumberOperators[Poly <: IPolynomial[Poly]](self: Long)(implicit ring: PolynomialRing[Poly, _]) {
    private val const: Poly = ring.valueOf(self)

    def +(poly: Poly): Poly = ring.add(const, poly)

    def -(poly: Poly): Poly = ring.subtract(const, poly)

    def *(poly: Poly): Poly = ring.multiply(const, poly)
  }

  implicit def operatorSupportU(number: Long)(implicit ring: PolynomialRing[UnivariatePolynomialZp64, Long])
  = NumberOperators[UnivariatePolynomialZp64](number)(ring)

  implicit def operatorSupportM(number: Long)(implicit ring: PolynomialRing[MultivariatePolynomialZp64, Long])
  = NumberOperators[MultivariatePolynomialZp64](number)(ring)

  implicit def operatorSupport[Poly <: IPolynomial[Poly]](number: Int)(implicit ring: PolynomialRing[Poly, _])
  = NumberOperators[Poly](number)(ring)

  /* ============================================= as coefficient ring ============================================= */

  implicit def asCoefficientDomainElementU[E](l: Long)(implicit ring: UnivariateRing[E]): E
  = ring.coefficientDomain.valueOf(l)

  implicit def asCoefficientDomainElementM[E](l: Long)(implicit ring: MultivariateRing[E]): E
  = ring.coefficientDomain.valueOf(l)

  implicit def asCoefficientDomainElementUi[E](l: Int)(implicit ring: UnivariateRing[E]): E
  = ring.coefficientDomain.valueOf(l)

  implicit def asCoefficientDomainElementMi[E](l: Int)(implicit ring: MultivariateRing[E]): E
  = ring.coefficientDomain.valueOf(l)

  implicit def asBigInteger(l: Int): Integer = BigInteger.valueOf(l)

  implicit def asBigInteger(l: Long): Integer = BigInteger.valueOf(l)

  def parse[Poly <: IPolynomial[Poly], E](string: String)(implicit ring: PolynomialRing[Poly, E])
  = ring.parse(string)

  def show[Poly <: IPolynomial[Poly], E](poly: Poly)(implicit ring: PolynomialRing[Poly, E])
  = scala.Predef.println(poly.toString(ring.variables))

  def println(arg: Any)(implicit ring: Ring[T] forSome {type T} = null): Unit = {
    scala.Predef.println(if (ring == null) arg.toString else ring.show(arg))
  }
}
