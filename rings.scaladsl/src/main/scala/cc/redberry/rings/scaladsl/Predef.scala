package cc.redberry.rings.scaladsl

import cc.redberry.rings
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.FiniteField
import cc.redberry.rings.{IntegersZp64, Rational, poly}

/**
  *
  * @since 1.0
  */
private[scaladsl] trait Predef {

  implicit def asBigInteger(v: BigInt): BigInteger = new BigInteger(v.bigInteger)

  implicit def asBigInteger(v: java.math.BigInteger): BigInteger = new BigInteger(v)

  implicit def asBigInteger(v: Int): BigInteger = BigInteger.valueOf(v)

  implicit def asBigInteger(v: Long): BigInteger = BigInteger.valueOf(v)

  implicit def asRingElement[E](v: Int)(implicit ring: Ring[E]) = ring.valueOf(v)

  implicit def asRingElement[E](v: Long)(implicit ring: Ring[E]) = ring.valueOf(v)

  /**
    * Delegate [[PolynomialRing]] methods for [[PolynomialRing]]
    */
  implicit def ringMethods[Poly <: IPolynomial[Poly], E](ring: PolynomialRing[Poly, E]): poly.PolynomialRing[Poly] = ring.theRing

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
  implicit def asRing[E](ring: rings.Ring[E]) = new Ring[E](ring)

  /**
    * Implicitly convert [[IntegersZp64]] to [[Ring]]
    */
  implicit def asRing(ring: IntegersZp64): Ring[Long] = new RingZp64(ring)

  /**
    * Ring of integers (Z)
    */
  val Z: Ring[BigInteger] = rings.Rings.Z

  /**
    * Field of rationals (Q)
    */
  val Q: Ring[Rational[BigInteger]] = rings.Rings.Q

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
             (implicit ring: PolynomialRing[MultivariatePolynomialZp64, Long]) = {
      val exps: Array[Int] = new Array[Int](exponents.length)
      exponents.foreach(e => exps(ring.index(e._1)) = e._2)
      new MonomialZp64(exps, ring.cfValue(cf))
    }
  }

  object Monomial {
    def apply[E](cf: E, exponents: Int*) = new Monomial[E](exponents.toArray, cf)

    def apply[E](cf: E, exponents: (String, Int)*)
             (implicit ring: PolynomialRing[MultivariatePolynomial[E], E]) = {
      val exps: Array[Int] = new Array[Int](exponents.length)
      exponents.foreach(e => exps(ring.index(e._1)) = e._2)
      new Monomial[E](exps, cf)
    }
  }
}
