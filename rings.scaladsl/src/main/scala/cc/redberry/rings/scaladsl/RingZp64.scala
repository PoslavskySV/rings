package cc.redberry.rings.scaladsl

import cc.redberry.rings
import cc.redberry.rings.IntegersZp64
import cc.redberry.rings.bigint.BigInteger

/**
  *
  * @since 1.0
  */
@specialized(Long)
private[scaladsl] final class RingZp64(val theRing: IntegersZp64) extends rings.Ring[Long] {
  override def isField = true

  override def isEuclideanRing = true

  override def cardinality() = BigInteger valueOf theRing.modulus

  override def characteristic() = BigInteger valueOf theRing.perfectPowerBase()

  override def isPerfectPower = theRing.isPerfectPower

  override def perfectPowerBase() = BigInteger valueOf theRing.perfectPowerBase()

  override def perfectPowerExponent() = BigInteger valueOf theRing.perfectPowerExponent()

  override def add(a: Long, b: Long) = theRing.add(a, b)

  override def subtract(a: Long, b: Long) = theRing.subtract(a, b)

  override def multiply(a: Long, b: Long) = theRing.multiply(a, b)

  override def negate(element: Long) = theRing.negate(element)

  override def copy(element: Long) = element

  override def divideAndRemainder(dividend: Long, divider: Long)
  = Array(theRing.divide(dividend, divider), 0L).asInstanceOf

  override def reciprocal(element: Long) = theRing.reciprocal(element)

  override def gcd(a: Long, b: Long) = a

  override def getZero = 0L

  override def getOne = 1L

  override def isZero(element: Long) = element == 0L

  override def isOne(element: Long) = element == 1L

  override def isUnit(element: Long) = true

  override def valueOfBigInteger(`val`: BigInteger) =
    if (`val`.isLong)
      valueOf(`val`.longValue())
    else
      `val`.mod(BigInteger valueOf theRing.modulus).longValueExact()

  override def valueOf(`val`: Long) = theRing.modulus(`val`)

  override def iterator() = Range(0, theRing.modulus.toInt).iterator.map(_.toLong).asInstanceOf

  override def compare(o1: Long, o2: Long) = java.lang.Long.compare(o1, o2)
}
