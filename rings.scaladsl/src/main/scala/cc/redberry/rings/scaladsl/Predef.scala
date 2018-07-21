package cc.redberry.rings.scaladsl

import cc.redberry.rings
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.FiniteField
import cc.redberry.rings.{FactorDecomposition, IntegersZp64, poly}
import org.apache.commons.math3.random.{AbstractRandomGenerator, RandomGenerator}

import scala.collection.JavaConverters.collectionAsScalaIterableConverter
import scala.language.implicitConversions
import scala.util.Random

/**
  *
  * @since 1.0
  */
private[scaladsl] trait Predef {

  implicit def asBigInteger(v: BigInt): IntZ = new BigInteger(v.bigInteger)

  implicit def asBigInteger(v: java.math.BigInteger): IntZ = new BigInteger(v)

  implicit def asBigInteger(v: Int): IntZ = BigInteger.valueOf(v)

  implicit def asBigInteger(v: Long): IntZ = BigInteger.valueOf(v)

  implicit def asRingElement[E](v: Int)(implicit ring: Ring[E]): E = ring.valueOf(v)

  implicit def asRingElement[E](v: Long)(implicit ring: Ring[E]) = ring.valueOf(v)

  implicit def asRandomGenerator(rnd: Random): RandomGenerator = new randomGenerator(rnd)

  private class randomGenerator(val random: Random) extends AbstractRandomGenerator {
    override def nextDouble(): Double = random.nextDouble()

    override def setSeed(seed: Long): Unit = random.setSeed(seed)
  }
  /**
    * Delegate [[IPolynomialRing]] methods for [[IPolynomialRing]]
    */
  implicit def ringMethods[Poly <: IPolynomial[Poly], E](ring: IPolynomialRing[Poly, E]): rings.poly.IPolynomialRing[Poly] = ring.theRing

  /**
    * Delegate [[IPolynomialRing]] methods for [[IPolynomialRing]]
    */
  implicit def ringMethods[E](ring: UnivariateRing[E]): poly.UnivariateRing[UnivariatePolynomial[E]]
  = ring.theRing.asInstanceOf[poly.UnivariateRing[UnivariatePolynomial[E]]]

  /**
    * Delegate [[IPolynomialRing]] methods for [[IPolynomialRing]]
    */
  implicit def ringMethods(ring: UnivariateRingZp64): poly.UnivariateRing[UnivariatePolynomialZp64]
  = ring.theRing.asInstanceOf[poly.UnivariateRing[UnivariatePolynomialZp64]]

  /**
    * Delegate [[IPolynomialRing]] methods for [[IPolynomialRing]]
    */
  implicit def ringMethods[E](ring: MultivariateRing[E]): poly.MultivariateRing[MultivariatePolynomial[E]]
  = ring.theRing.asInstanceOf[poly.MultivariateRing[MultivariatePolynomial[E]]]

  /**
    * Delegate [[IPolynomialRing]] methods for [[IPolynomialRing]]
    */
  implicit def ringMethods(ring: MultivariateRingZp64): poly.MultivariateRing[MultivariatePolynomialZp64]
  = ring.theRing.asInstanceOf[poly.MultivariateRing[MultivariatePolynomialZp64]]

  /**
    * Delegate [[FiniteField]] methods fo [[GaloisField64]]
    */
  implicit def ringMethods(ring: GaloisField64): FiniteField[UnivariatePolynomialZp64] = ring.theRing

  /**
    * Delegate [[FiniteField]] methods fo [[GaloisField64]]
    */
  implicit def ringMethods[E](ring: GaloisField[E]): FiniteField[UnivariatePolynomial[E]] = ring.theRing

  /**
    * Delegate [[rings.Ring]] methods for [[Ring]]
    */
  implicit def ringMethods[E](ring: Ring[E]): rings.Ring[E] = ring.theRing

  /**
    * Implicitly convert [[rings.Ring]] to [[Ring]]
    */
  implicit def asRing[E](ring: rings.Ring[E]): Ring[E] = new Ring[E](ring)

  /**
    * Implicitly convert [[IntegersZp64]] to [[Ring]]
    */
  implicit def asRing(ring: IntegersZp64): Ring[Long] = new RingZp64(ring)

  /**
    * Delegate [[Ideal]] methods for [[Ideal]]
    */
  implicit def idealMethods[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](ideal: Ideal[Term, Poly, E]):
  rings.poly.multivar.Ideal[Term, Poly] = ideal.theIdeal

  /**
    * Implicitly convert [[rings.poly.multivar.Ideal]] to [[Ideal]]
    */
  implicit def asIdeal[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], E](ideal: poly.multivar.Ideal[Term, Poly]):
  Ideal[Term, Poly, E] = new Ideal[Term, Poly, E](MultivariateRing[Term, Poly, E](ideal.getBasisGenerator(0)), ideal)

  /**
    * Ring of integers (Z)
    */
  val Z: Ring[BigInteger] = rings.Rings.Z

  /**
    * Field of rationals (Q)
    */
  val Q: Ring[rings.Rational[BigInteger]] = rings.Rings.Q

  /**
    * Field of integers modulo `modulus`
    *
    * @param modulus the modulus
    */
  def Zp64(modulus: Long): rings.IntegersZp64 = rings.Rings.Zp64(modulus)

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
    * Field of integers modulo `modulus`
    *
    * @param modulus the modulus
    */
  def Zp(modulus: BigInt): Ring[BigInteger] = Zp(new BigInteger(modulus.bigInteger))

  /**
    * Gaussian numbers for a given ring (that is ring adjoined with imaginary unit)
    */
  def GaussianNumbers[E](ring: Ring[E], imaginaryUnit: String = "i")
  : AlgebraicNumberField[E] = AlgebraicNumberField(rings.Rings.GaussianNumbers(ring), imaginaryUnit)

  /**
    * Ring of Gaussian integers (integer complex numbers).
    */
  lazy val GaussianIntegers: AlgebraicNumberField[IntZ] = GaussianIntegers("i")

  /**
    * Ring of Gaussian integers (integer complex numbers).
    */
  def GaussianIntegers(imaginaryUnit: String = "i")
  : AlgebraicNumberField[IntZ] = AlgebraicNumberField(rings.Rings.GaussianIntegers, imaginaryUnit)

  /**
    * Field of Gaussian rationals (rational complex numbers).
    */
  lazy val GaussianRationals: AlgebraicNumberField[Rational[IntZ]] = GaussianRationals("i")

  /**
    * Field of Gaussian rationals (rational complex numbers).
    */
  def GaussianRationals(imaginaryUnit: String)
  : AlgebraicNumberField[Rational[IntZ]] = AlgebraicNumberField(rings.Rings.GaussianRationals, imaginaryUnit)

  /**
    * Splitting field of a given polynomial.
    */
  def SplittingField[
  Term <: AMonomial[Term],
  mPoly <: AMultivariatePolynomial[Term, mPoly],
  sPoly <: IUnivariatePolynomial[sPoly],
  E](poly: sPoly, variables: Array[String])
  : MultipleFieldExtension[Term, mPoly, sPoly, E]
  = MultipleFieldExtension(rings.Rings.SplittingField(poly), variables)

  /**
    * Splitting field of a given polynomial.
    */
  def SplittingField[E](poly: UnivariatePolynomial[E], variables: Array[String])
  : MultipleFieldExtension[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E], E]
  = MultipleFieldExtension(rings.Rings.SplittingField[Monomial[E], MultivariatePolynomial[E], UnivariatePolynomial[E]](poly), variables)

  /**
    * Splitting field of a given polynomial.
    */
  def SplittingFieldZp64(poly: UnivariatePolynomialZp64, variables: Array[String])
  : MultipleFieldExtension[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64, Long]
  = MultipleFieldExtension(rings.Rings.SplittingField[MonomialZp64, MultivariatePolynomialZp64, UnivariatePolynomialZp64](poly), variables)

  object UnivariatePolynomial {
    def apply[I, E](cfs: I*)(implicit ring: Ring[E]): UnivariatePolynomial[E] = cfs.headOption match {
      case None => rings.poly.univar.UnivariatePolynomial.zero[E](ring)
      case Some(head) => head match {
        case _: Int =>
          val longs: Array[Long] = cfs.map(_.asInstanceOf[Int].toLong).toArray
          rings.poly.univar.UnivariatePolynomial.create[E](ring, ring.valueOf(longs): _*)
        case h: E =>
          val es: Array[E] = ring.createArray(cfs.length).asInstanceOf[Array[E]]
          cfs.map(e => ring.valueOf(e.asInstanceOf[E])).copyToArray(es)
          rings.poly.univar.UnivariatePolynomial.create[E](ring, es: _*)
        case _ => throw new IllegalArgumentException(s"unknown type $head.getC")
      }
    }
  }

  object UnivariatePolynomialZp64 {
    def apply(cfs: Long*)(implicit ring: RingZp64): UnivariatePolynomialZp64
    = rings.poly.univar.UnivariatePolynomialZp64.create(ring.theRing, cfs.toArray)
  }

  object MonomialZp64 {
    def apply(cf: Int, exponents: Int*) = new MonomialZp64(exponents.toArray, cf.toLong)

    def apply(cf: Int, exponents: (String, Int)*)
             (implicit ring: IPolynomialRing[MultivariatePolynomialZp64, Long]) = {
      val exps: Array[Int] = new Array[Int](exponents.length)
      exponents.foreach(e => exps(ring.index(e._1)) = e._2)
      new MonomialZp64(exps, ring.cfValue(cf))
    }
  }

  object Monomial {
    def apply[E](cf: E, exponents: Int*) = new Monomial[E](exponents.toArray, cf)

    def apply[E](cf: E, exponents: (String, Int)*)
                (implicit ring: IPolynomialRing[MultivariatePolynomial[E], E]) = {
      val exps: Array[Int] = new Array[Int](exponents.length)
      exponents.foreach(e => exps(ring.index(e._1)) = e._2)
      new Monomial[E](exps, cf)
    }
  }

  object Rational {
    def apply[E](num: E, den: E)(implicit ring: Ring[E]): rings.Rational[E] = new rings.Rational[E](ring, num, den)

    def apply[E](num: E)(implicit ring: Ring[E]): rings.Rational[E] = new rings.Rational[E](ring, num)

    @inline
    def apply[E](num: Int, den: E)(implicit ring: Ring[E]): rings.Rational[E] = apply(ring(num), den)(ring)

    @inline
    def apply[E](num: E, den: Int)(implicit ring: Ring[E]): rings.Rational[E] = apply(num, ring(den))(ring)

    @inline
    def apply[E](num: Int, den: Int)(implicit ring: Ring[E]): rings.Rational[E] = apply(ring(num), ring(den))(ring)

    @inline
    def apply[E](num: Int)(implicit ring: Ring[E]): rings.Rational[E] = apply(ring(num))(ring)
  }

  final class RichArrayTuple[Poly](arr: Array[Poly]) {
    def tuple2: (Poly, Poly) = (arr(0), arr(1))

    def tuple3: (Poly, Poly, Poly) = (arr(0), arr(1), arr(2))

    def tuple4: (Poly, Poly, Poly, Poly) = (arr(0), arr(1), arr(2), arr(3))
  }

  implicit def arrayToTuple[Poly](arr: Array[Poly]): RichArrayTuple[Poly] = new RichArrayTuple(arr)

  implicit def factors2Seq[E](factors: FactorDecomposition[E]): Seq[(E, Int)] =
    (factors.factors.asScala zip factors.exponents.toArray().toSeq).toSeq
}
