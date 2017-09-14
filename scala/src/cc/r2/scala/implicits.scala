package cc.r2.scala

import cc.r2.core.number.BigInteger
import cc.r2.core.poly.multivar._
import cc.r2.core.poly.univar.{IUnivariatePolynomial, UnivariatePolynomial, UnivariatePolynomialZp64}
import cc.r2.core.poly.{IPolynomial, PolynomialMethods}

import scala.reflect.ClassTag

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
object implicits {
  type Integer = BigInteger

  import PolynomialRings._

  private[implicits] sealed class IPolynomialOperators[Poly <: IPolynomial[Poly]]
  (self: Poly)(implicit ring: PolynomialRing[Poly]) {
    protected[implicits] def constant(value: Int): Poly = self.createConstant(value)

    protected[implicits] def constant(value: Long): Poly = self.createConstant(value)

    def +(other: Poly): Poly = ring.add(self, other)

    def +(other: Long): Poly = this.+(constant(other))

    def -(other: Poly): Poly = ring.subtract(self, other)

    def -(other: Long): Poly = this.-(constant(other))

    def *(other: Poly): Poly = ring.multiply(self, other)

    def *(other: Long): Poly = this.*(constant(other))

    def **(exponent: Int): Poly = PolynomialMethods.polyPow(self, exponent, true)

    def ^(exponent: Int): Poly = this.**(exponent)

    def unary_- : Poly = ring.negate(self)

    def unary_+ : Poly = self.copy()

    def ++ : Poly = ring.increment(self)

    def -- : Poly = ring.decrement(self)

    def <(other: Poly): Boolean = self.compareTo(other) < 0

    def <=(other: Poly): Boolean = self.compareTo(other) <= 0

    def >(other: Poly): Boolean = self.compareTo(other) > 0

    def >=(other: Poly): Boolean = self.compareTo(other) >= 0

    def /%(other: Poly): (Poly, Poly) = {
      val qd = ring.divideAndRemainder(self, other)
      (qd(0), qd(1))
    }

    def /(other: Poly): Poly = ring.divideExact(self, other)

    def /(other: Long): Poly = this./(constant(other))

    def ===(other: Int): Boolean = self == constant(other)

    def ===(other: Poly): Boolean = self == other
  }

  private[implicits] sealed class IUnivariateOperators[Poly <: IUnivariatePolynomial[Poly]](self: Poly)(implicit ring: PolynomialRing[Poly])
    extends IPolynomialOperators[Poly](self)(ring) {

    def <<(offset: Int): Poly = ring.valueOf(self.copy().shiftLeft(offset))

    def >>(offset: Int): Poly = ring.valueOf(self.copy().shiftRight(offset))

    def apply(from: Int, to: Int): Poly = ring.valueOf(self.copy().getRange(from, to))

    def apply(poly: Poly): Poly = ring.valueOf(self.composition(poly))
  }

  final case class UnivariateZp64Operators
  (self: UnivariatePolynomialZp64)(implicit ring: PolynomialRing[UnivariatePolynomialZp64])
    extends IUnivariateOperators[UnivariatePolynomialZp64](self)(ring) {
    def apply(index: Int): Long = self.get(index)

    def eval(point: Long): Long = self.evaluate(point)
  }

  final case class UnivariateOperators[ElementType]
  (self: UnivariatePolynomial[ElementType])(implicit ring: PolynomialRing[UnivariatePolynomial[ElementType]])
    extends IUnivariateOperators[UnivariatePolynomial[ElementType]](self)(ring) {
    def apply(index: Int): ElementType = self.get(index)

    def eval(point: ElementType): ElementType = self.evaluate(point)
  }

  private[implicits] sealed abstract class IMultivariateOperators[Term <: DegreeVector[Term], Poly <: AMultivariatePolynomial[Term, Poly]]
  (self: Poly)(implicit ring: PolynomialRing[Poly])
    extends IPolynomialOperators[Poly](self)(ring) {
    def +(other: Term): Poly = ring.pDomain.valueOf(self.copy().add(other))

    def -(other: Term): Poly = ring.pDomain.valueOf(self.copy().subtract(other))

    def *(other: Term): Poly = ring.pDomain.valueOf(self.copy().multiply(other))

    def /(other: Term): Poly = {
      val r = self.copy().divideOrNull(other)
      if (r == null)
        throw new ArithmeticException(s"not divisible: $self / $other")
      else
        ring.pDomain.valueOf(r)
    }

    def swapVariables(i: Int, j: Int): Poly

    def swapVariables(tup: (String, String)): Poly = swapVariables(ring.index(tup._1), ring.index(tup._2))
  }

  final case class MultivariateZp64Operators
  (self: MultivariatePolynomialZp64)(implicit ring: PolynomialRing[MultivariatePolynomialZp64])
    extends IMultivariateOperators[MonomialZp64, MultivariatePolynomialZp64](self)(ring) {

    override def swapVariables(i: Int, j: Int) = AMultivariatePolynomial.swapVariables[MonomialZp64, MultivariatePolynomialZp64](self, i, j)

    def eval(variable: String, value: Long): MultivariatePolynomialZp64 = self.evaluate(ring.index(variable), value)

    def eval(tup: (String, Long)): MultivariatePolynomialZp64 = eval(tup._1, tup._2)

    def eval(tups: (String, Long)*): MultivariatePolynomialZp64 = self.evaluate(tups.map(t => ring.index(t._1)).toArray, tups.map(_._2).toArray)
  }

  final case class MultivariateOperators[E <: AnyRef]
  (self: MultivariatePolynomial[E])(implicit ring: PolynomialRing[MultivariatePolynomial[E]])
    extends IMultivariateOperators[Monomial[E], MultivariatePolynomial[E]](self)(ring) {

    override def swapVariables(i: Int, j: Int) = AMultivariatePolynomial.swapVariables[Monomial[E], MultivariatePolynomial[E]](self, i, j)

    def eval(variable: String, value: Long): MultivariatePolynomial[E] = self.evaluate(ring.index(variable), value)

    def eval(tup: (String, Long)): MultivariatePolynomial[E] = eval(tup._1, tup._2)

    private implicit val eClassTag: ClassTag[E] = ClassTag(self.cc().getClass)

    def eval(tups: (String, E)*): MultivariatePolynomial[E] = self.evaluate1(tups.map(t => ring.index(t._1)).toArray[Int], tups.map(_._2).toArray[E])
  }

  final case class StringSwap(self: String) {
    def <>(oth: String) = (self, oth)

    def ->(oth: Long) = (self, oth)

    def ->[E](oth: E) = (self, oth)
  }

  implicit def stringWithSwap(self: String) = StringSwap(self)

  implicit def operatorSupport
  (poly: UnivariatePolynomialZp64)
  (implicit ring: PolynomialRing[UnivariatePolynomialZp64]) =
    UnivariateZp64Operators(poly)(ring)

  implicit def operatorSupport[ElementType]
  (poly: UnivariatePolynomial[ElementType])
  (implicit ring: PolynomialRing[UnivariatePolynomial[ElementType]]) =
    UnivariateOperators[ElementType](poly)(ring)

  implicit def operatorSupport
  (poly: MultivariatePolynomialZp64)
  (implicit ring: PolynomialRing[MultivariatePolynomialZp64]) =
    MultivariateZp64Operators(poly)(ring)

  implicit def operatorSupport[ElementType <: AnyRef]
  (poly: MultivariatePolynomial[ElementType])
  (implicit ring: PolynomialRing[MultivariatePolynomial[ElementType]]) =
    MultivariateOperators(poly)(ring)

  final case class NumberOperators[Poly <: IPolynomial[Poly]](self: Long)(implicit ring: PolynomialRing[Poly]) {
    private implicit def operatorSupport(poly: Poly) = new IPolynomialOperators(poly)

    def +(poly: Poly): Poly = poly.+(self)

    def -(poly: Poly): Poly = -(poly - self)

    def *(poly: Poly): Poly = poly * self
  }

  implicit def operatorSupport[Poly <: IPolynomial[Poly]]
  (number: Long)
  (implicit ring: PolynomialRing[Poly]) =
    NumberOperators(number)(ring)


  implicit def asBigInteger(l: Long) = BigInteger.valueOf(l)

  def parse[PolynomialType <: IPolynomial[PolynomialType]](string: String)(implicit ring: PolynomialRing[PolynomialType])
  = ring.parse(string)

  def show[PolynomialType <: IPolynomial[PolynomialType]](poly: PolynomialType)(implicit ring: PolynomialRing[PolynomialType])
  = scala.Predef.println(poly.toString(ring.variables: _*))

  def println(any: Any)(implicit ring: PolynomialRing[_]) {
    any match {
      case poly: IPolynomial[_] =>
        scala.Predef.println(poly.toString(ring.variables: _*))
      case _ => scala.Predef.println(any)
    }
  }
}
