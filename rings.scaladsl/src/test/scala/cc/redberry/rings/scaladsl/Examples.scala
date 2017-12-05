package cc.redberry.rings.scaladsl

import java.util

import cc.redberry.rings
import cc.redberry.rings.Rational
import cc.redberry.rings.linear.LinearSolver
import cc.redberry.rings.poly.PolynomialMethods
import cc.redberry.rings.poly.multivar.MonomialOrder.LEX
import cc.redberry.rings.poly.multivar.{MonomialOrder, MultivariateInterpolation, MultivariatePolynomial}
import cc.redberry.rings.poly.univar.{IrreduciblePolynomials, UnivariateGCD, UnivariatePolynomialArithmetic}
import cc.redberry.rings.primes.SmallPrimes
import org.apache.commons.math3.random.Well1024a
import org.junit.{Ignore, Test}

import scala.collection.mutable
import scala.language.{existentials, postfixOps}
import scala.util.Random

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
class Examples {
  @Test
  def test1: Unit = {

    import syntax._
    // when parsing "x" will be considered as the "first variable"
    // and "y" as "the second" => in the result the particular
    // names "x" and "y" are erased
    val poly1 = MultivariatePolynomial.parse("x^2 + x*y", "x", "y")
    // parse the same polynomial but using "a" and "b" instead of "x" and "y"
    val poly2 = MultivariatePolynomial.parse("a^2 + a*b", "a", "b")
    // polynomials are equal (no matter which variable names were used when parsing)
    assert(poly1 == poly2)
    // degree in the first variable
    assert(poly1.degree(0) == 2)
    // degree in the second variable
    assert(poly1.degree(1) == 1)

    // this poly differs from poly2 since now "a" is "the second"
    // variable and "b" is "the first"
    val poly3 = MultivariatePolynomial.parse("a^2 + a*b", "b", "a")
    assert(poly3 != poly2)
    // swap the first and the second variables and the result is equal to poly2
    assert(poly3.swapVariables(0, 1) == poly2)


    // the default toString() will use the default
    // variables "a", "b", "c"  and so on (alphabetical)
    // the result will be "a*b + a^2"
    println(poly1)
    // specify which variable names use for printing
    // the result will be "x*y + x^2"
    println(poly1.toString(Array("x", "y")))
    // the result will be "y*x + y^2"
    println(poly1.toString(Array("y", "x")))
  }

  @Test
  def test2: Unit = {
    // "x" is the first variable "y" is the second
    val ring = MultivariateRing(Z, Array("x", "y"))
    // parse polynomial
    val poly = ring("x^2 + x*y")
    // the result will be "x*y + x^2"
    println(ring show poly)
  }

  @Test
  def test3: Unit = {

    val ring = UnivariateRing(Z, "x")
    // parse polynomials
    val (p1, p2, p3) = ring("x", "x^2", "x^3")

    // this WILL modify poly1
    p1.add(p2)
    // this will NOT modify poly2
    p2.copy().add(p3)
    // this will NOT modify poly2
    ring.add(p2, p3)
  }

  @Test
  def test4: Unit = {
    import syntax._

    val ring = Q

    val a = ring("2/3")
    val b = ring("3/4")

    println(a + b)
  }

  @Test
  def test5: Unit = {
    import syntax._

    implicit val ring = UnivariateRing(Zp(3), "x")
    val (a, b) = ring("1 + 2*x^2", "1 - x")

    // add two elements
    val add = a + b
    // subtract two elements
    val sub = a - b
    // multiply two elements
    val mul = a * b
    // divide two elements
    val div = a / b
    // divide with remainder
    val divRem = a /% b
    // add one
    val inc = a ++
    // subtract one
    val dec = a --
    // negate element
    val neg = -a


  }

  @Test
  def test6: Unit = {
    import syntax._

    val a: IntZ = 10
    val b: IntZ = 11

    // compiles to a.add(b) (integer addition)
    assert(a + b === 21)


    implicit val ring = Zp(13)
    // compiles to ring.add(a, b) (addition mod 13)
    assert(a + b === 8)
  }

  @Test
  def test7: Unit = {
    import syntax._

    implicit val ring = UnivariateRingZp64(17, "x")
    // some random divider
    val divider = ring.randomElement()
    // some random dividend
    val dividend = 1 + 2 * divider + 3 * divider.pow(2)

    // quotient and remainder using built-in methods
    val (divPlain, remPlain) = dividend /% divider

    // precomputed Newton inverses, need to calculate it only once
    implicit val invMod: PrecomputedInverse[UnivariatePolynomialZp64] = divider.precomputedInverses
    // quotient and remainder computed using fast
    // algorithm with precomputed Newton inverses
    val (divFast, remFast) = dividend /%% divider

    // result is the same
    assert((divPlain, remPlain) == (divFast, remFast))
  }

  @Test
  def test8: Unit = {
    import rings.poly.univar.UnivariateGCD._
    import syntax._

    // Polynomials over field
    val ringZp = UnivariateRingZp64(17, "x")
    val a = ringZp("1 + 3*x + 2*x^2")
    val b = ringZp("1 - x^2")
    // Euclid and Half-GCD algorithms for polynomials over field
    assert(EuclidGCD(a, b) == HalfGCD(a, b))
    // Extended Euclidean algorithm
    val (gcd, s, t) = ExtendedEuclidGCD(a, b) match {
      case Array(gcd, s, t) => (gcd, s, t)
    }
    assert(a * s + b * t == gcd)
    // Extended Half-GCD algorithm
    val (gcd1, s1, t1) = ExtendedHalfGCD(a, b) match {
      case Array(gcd, s, t) => (gcd, s, t)
    }
    assert((gcd1, s1, t1) == (gcd, s, t))


    // Polynomials over Z
    val ringZ = UnivariateRing(Z, "x")
    val aZ = ringZ("1 + 3*x + 2*x^2")
    val bZ = ringZ("1 - x^2")
    // GCD for polynomials over Z
    assert(ModularGCD(aZ, bZ) == ringZ("1 + x"))

    // Bivariate polynomials represented as Z[y][x]
    val ringXY = UnivariateRing(UnivariateRing(Z, "y"), "x")
    val aXY = ringXY("(1 + y) + (1 + y^2)*x + (y - y^2)*x^2")
    val bXY = ringXY("(3 + y) + (3 + 2*y + y^2)*x + (3*y - y^2)*x^2")
    // Subresultant sequence
    val subResultants = SubresultantRemainders(aXY, bXY)
    // The GCD
    val gcdXY = subResultants.gcd.primitivePart
    assert(aXY % gcdXY === 0 && bXY % gcdXY === 0)
  }


  @Test
  def test9: Unit = {
    import rings.poly.PolynomialMethods._
    import rings.poly.univar.UnivariateSquareFreeFactorization._
    import rings.scaladsl.syntax._

    // ring GF(13^5)[x] (coefficient domain is finite field)
    val ringF = UnivariateRing(GF(13, 5, "z"), "x")
    // some random polynomial composed from some factors
    val polyF = ringF.randomElement() * ringF.randomElement() * ringF.randomElement().pow(10)
    // perform square-free factorization
    println(ringF show SquareFreeFactorization(polyF))
    // perform complete factorization
    println(ringF show Factor(polyF))


    // ring Q[x]
    val ringQ = UnivariateRing(Q, "x")
    // some random polynomial composed from some factors
    val polyQ = ringQ.randomElement() * ringQ.randomElement() * ringQ.randomElement().pow(10)
    // perform square-free factorization
    println(ringQ show SquareFreeFactorization(polyQ))
    // perform complete factorization
    println(ringQ show Factor(polyQ))
  }


  @Test
  def test10: Unit = {
    import rings.poly.univar.IrreduciblePolynomials._
    val random = new Random()

    // random irreducible polynomial in Z/2[x] of degree 10
    val poly1: UnivariatePolynomialZp64 = randomIrreduciblePolynomial(2, 10, random)
    assert(poly1.degree() == 10)
    assert(irreducibleQ(poly1))

    // random irreducible polynomial in Z/2[x] of degree 10
    val poly2: UnivariatePolynomial[IntZ] = randomIrreduciblePolynomial(Zp(2).theRing, 10, random)
    assert(poly2.degree() == 10)
    assert(irreducibleQ(poly2))

    // random irreducible polynomial in GF(11^15)[x] of degree 10
    val poly3: UnivariatePolynomial[UnivariatePolynomialZp64] = randomIrreduciblePolynomial(GF(11, 15).theRing, 10, random)
    assert(poly3.degree() == 10)
    assert(irreducibleQ(poly3))

    // random irreducible polynomial in Z[x] of degree 10
    val poly4: UnivariatePolynomial[IntZ] = randomIrreduciblePolynomialOverZ(10, random)
    assert(poly4.degree() == 10)
    assert(irreducibleQ(poly4))
  }

  @Test
  def test11: Unit = {
    import rings.poly.univar.UnivariateInterpolation._

    // points
    val points = Array(1L, 2L, 3L, 12L)
    // values
    val values = Array(3L, 2L, 1L, 6L)

    val result = new InterpolationZp64(Zp64(17))
      .update(points, values)
      .getInterpolatingPolynomial

    assert(points.zipWithIndex.forall { case (point, i) => result.evaluate(point) == values(i) })

  }


  @Test
  def test12: Unit = {

    /*  Lagrange interpolation formula */
    def lagrange[Poly <: IUnivariatePolynomial[Poly], E](points: Seq[E], values: Seq[E])(implicit ring: IUnivariateRing[Poly, E]) = {
      import syntax._
      points.indices
        .foldLeft(ring getZero) { case (sum, i) =>
          sum + points.indices
            .filter(_ != i)
            .foldLeft(ring getConstant values(i)) { case (product, j) =>
              implicit val cfRing = ring.cfRing
              val E: E = points(i) - points(j)
              product * (ring.`x` - points(j)) / E
            }
        }
    }

    import rings.poly.univar.UnivariateInterpolation._
    import syntax._

    // coefficient ring GF(13, 5)
    implicit val cfRing = GF(13, 5, "z")
    val z = cfRing("z")
    // some points
    val points = Array(1 + z, 2 + z, 3 + z, 12 + z)
    // some values
    val values = Array(3 + z, 2 + z, 1 + z, 6 + z)

    // interpolate with Newton iterations
    val withNewton = new Interpolation(cfRing)
      .update(points, values)
      .getInterpolatingPolynomial
    // interpolate using Lagrange formula
    val withLagrange = lagrange(points, values)(UnivariateRing(cfRing, "x"))
    // results are the same
    assert(withNewton == withLagrange)
  }

  @Test
  def test13 = {
    /*  Lagrange interpolation formula */
    def lagrange[Poly <: IUnivariatePolynomial[Poly], E](points: Seq[E], values: Seq[E])(implicit ring: IUnivariateRing[Poly, E]) = {
      import syntax._
      points.indices
        .foldLeft(ring getZero) { case (sum, i) =>
          sum + points.indices
            .filter(_ != i)
            .foldLeft(ring getConstant values(i)) { case (product, j) =>
              implicit val cfRing = ring.cfRing
              val E: E = points(i) - points(j)
              product * (ring.`x` - points(j)) / E
            }
        }
    }
  }

  @Test
  def test14 = {

    /**
      * @tparam Poly type of polynomials
      */
    def genericFunc[Poly <: IPolynomial[Poly]](poly: Poly): Poly = {
      import syntax._
      poly.pow(2) * 3 + poly * 2 + 1
    }

    // univariate polynomials over Zp64
    val uRing = UnivariateRingZp64(17, "x")
    println(uRing show genericFunc(uRing("1 + 2*x + 3*x^2")))

    // multivariate polynomials over Z
    val mRing = MultivariateRing(Z, Array("x", "y", "z"))
    println(mRing show genericFunc(mRing("1 + x + y + z")))
  }

  @Test
  def test15 = {

    /**
      * @tparam Poly type of polynomials
      * @tparam E    type of polynomial coefficients
      */
    def genericFuncWithRing[Poly <: IPolynomial[Poly], E](poly: Poly)
                                                         (implicit ring: IPolynomialRing[Poly, E]): Poly = {
      import syntax._
      poly.pow(2) * 3 + poly * 2 + 1
    }

    // univariate polynomials over Zp64
    val uRing = UnivariateRingZp64(17, "x")
    println(uRing show genericFuncWithRing(uRing("1 + 2*x + 3*x^2"))(uRing))

    // multivariate polynomials over Z
    val mRing = MultivariateRing(Z, Array("x", "y", "z"))
    println(mRing show genericFuncWithRing(mRing("1 + x + y + z"))(mRing))
  }


  @Test
  def test16 = {
    /**
      * @tparam Poly type of polynomials
      * @tparam E    type of polynomial coefficients
      */
    def genericFunc[Poly <: IPolynomial[Poly], E](poly: Poly): Poly = {
      import syntax._
      poly.pow(2) * 3 + poly * 2 + 1
    }

    /**
      * @tparam Poly type of polynomials
      * @tparam E    type of polynomial coefficients
      */
    def genericFuncWithRing[Poly <: IPolynomial[Poly], E](poly: Poly)
                                                         (implicit ring: IPolynomialRing[Poly, E]): Poly = {
      import syntax._
      poly.pow(2) * 3 + poly * 2 + 1
    }

    import syntax._
    // GF(13^4)
    implicit val gf = GF(13, 4, "z")
    // some element of GF(13^4)
    val poly = gf("1 + z + z^2 + z^3 + z^4").pow(10)

    val noring = genericFunc(poly)
    println(noring)

    val withring = genericFuncWithRing(poly)
    println(withring)

    assert(noring != withring)
  }


  @Test
  def test17 = {
    /**
      * @tparam Poly type of univariate polynomials
      * @tparam E    type of polynomial coefficients
      */
    def genericFunc[Poly <: IUnivariatePolynomial[Poly], E](poly: Poly) = poly

    /**
      * @tparam Poly type of univariate polynomials
      * @tparam E    type of polynomial coefficients
      */
    def genericFuncWithRing[Poly <: IUnivariatePolynomial[Poly], E](poly: Poly)
                                                                   (implicit ring: IUnivariateRing[Poly, E]) = poly
  }

  @Test
  def test18: Unit = {
    /**
      * @tparam Monomial type of monomials
      * @tparam Poly     type of multivariate polynomials
      * @tparam E        type of polynomial coefficients
      */
    def genericFunc[
    Monomial <: DegreeVector[Monomial],
    Poly <: AMultivariatePolynomial[Monomial, Poly],
    E
    ](poly: Poly) = poly

    /**
      * @tparam Monomial type of monomials
      * @tparam Poly     type of multivariate polynomials
      * @tparam E        type of polynomial coefficients
      */
    def genericFuncWithRing[Monomial <: DegreeVector[Monomial], Poly <: AMultivariatePolynomial[Monomial, Poly], E](poly: Poly)
                                                                                                                   (implicit ring: IMultivariateRing[Monomial, Poly, E]) = poly


    implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))
    import ring.{CoefficientType, MonomialType, PolyType}

    val poly = ring.randomElement()

    // call generic func directly
    genericFunc[MonomialType, PolyType, CoefficientType](poly)
    genericFuncWithRing[MonomialType, PolyType, CoefficientType](poly)

    // define shortcuts
    val func = (p: ring.PolyType) =>
      genericFunc[MonomialType, PolyType, CoefficientType](p)

    val funcWithRing = (p: ring.PolyType) => genericFuncWithRing[MonomialType, PolyType, CoefficientType](p)(ring)

    func(poly)
    funcWithRing(poly)

    LEX
  }


  @Test
  def test19: Unit = {
    import syntax._
    val ring = MultivariateRing(Z, Array("x", "y"), MonomialOrder.GREVLEX)

    // poly in GREVLEX
    val poly = ring("x + x^2*y^2 + x*y")
    assert(poly.ordering == MonomialOrder.GREVLEX)

    // poly in LEX
    val poly2 = poly.setOrdering(LEX)
    assert(poly2.ordering == LEX)

    // poly in GREVLEX (ordering of lhs is used)
    val add = poly + poly2
    assert(add.ordering == MonomialOrder.GREVLEX)

    // poly in LEX (ordering of lhs is used)
    val add2 = poly2 + poly
    assert(add2.ordering == LEX)
  }

  @Test
  def test20: Unit = {
    import syntax._
    val ring = MultivariateRing(Z, Array("x", "y", "z"), LEX)

    val dividend = ring("x - x^2*y^2 + 2*x*y + 1 - z*y^2*x^2 + z").pow(3)

    val divider1 = ring("x + y")
    val divider2 = ring("x + z")
    val divider3 = ring("y + z")

    {
      val (quot1, quot2, rem) = dividend /%/% (divider1, divider2)
      assert(dividend == divider1 * quot1 + divider2 * quot2 + rem)
    }


    {
      val (quot1, quot2, quot3, rem) = dividend /%/% (divider1, divider2, divider3)
      assert(dividend == divider1 * quot1 + divider2 * quot2 + divider3 * quot3 + rem)
    }
  }

  @Test
  def test21: Unit = {
    import rings.poly.multivar.MultivariateGCD._
    import syntax._


    // some large finite field
    val modulus = SmallPrimes.nextPrime(1 << 15)
    val ring = MultivariateRingZp64(modulus, Array("x", "y", "z"))

    val a = ring("x^2 - x*y + z^5")
    val b = ring("x^2 + x*y^7 + x*y*z^2")

    val gcd = ring("x + y + z")
    val poly1 = a * gcd
    val poly2 = b * gcd

    // EZGCD in finite field
    val ez = EZGCD(poly1, poly2)
    assert(ez == gcd)

    // EEZGCD in finite field
    val eez = EEZGCD[ring.MonomialType, ring.PolyType](poly1, poly2)
    assert(eez == gcd)

    // ZippelGCD in finite field
    val zippel = ZippelGCD(poly1, poly2)
    assert(zippel == gcd)

    // some very small finite field (Z/2)
    val z2 = Zp64(2)
    val z2GCD = gcd.setRing(z2)
    val z2Poly1 = a.setRing(z2) * z2GCD
    val z2Poly2 = b.setRing(z2) * z2GCD

    // Kaltofen’s & Monagan’s generic modular GCD
    val modGF = KaltofenMonaganSparseModularGCDInGF(z2Poly1, z2Poly2)
    assert(modGF == z2GCD)

    // Z
    val zGCD = gcd.setRing[IntZ](Z)
    val zPoly1 = a.setRing[IntZ](Z) * zGCD
    val zPoly2 = b.setRing[IntZ](Z) * zGCD

    // Modular GCD in Z with sparse interpolation
    val mod = ZippelGCDInZ(zPoly1, zPoly2)
    assert(mod == zGCD)
  }


  @Test
  def test22: Unit = {
    import rings.poly.multivar.MultivariateGCD._
    import syntax._

    val ring = MultivariateRing(Z, Array("x", "y", "z"))
    val (rndDegree, rndSize) = (5, 5)

    // some random gcd
    val gcd = ring.randomElement(rndDegree, rndSize)
    // array of random polynomials which have gcd
    val polys = (0 until 10).map(_ => ring.randomElement(rndDegree, rndSize) * gcd)

    // fast algorithm for array of polynomials will be used
    val fastGCD = PolynomialGCD(polys: _*)
    // slow step-by-step gcd calculation
    val slowGCD = polys.foldLeft(ring.getZero)((p1, p2) => PolynomialGCD(p1, p2))
    // result the same (to within a sign)
    assert(fastGCD ~== slowGCD)

  }

  @Ignore
  @Test
  def test23: Unit = {
    import rings.poly.PolynomialMethods._
    import syntax._

    // ring GF(13^5)[x, y, z] (coefficient domain is finite field)
    val ringF = MultivariateRing(GF(13, 5), Array("x", "y", "z"))

    // generate random poly of degree 5 and size 5
    def randomPolyF = ringF.randomElement(5, 5) + 1

    // some random polynomial composed from some factors
    val polyF = randomPolyF * randomPolyF * randomPolyF.pow(2)
    // perform square-free factorization
    println(ringF show FactorSquareFree(polyF))
    // perform complete factorization
    println(ringF show Factor(polyF))


    // ring Q[x, y, z]
    val ringQ = MultivariateRing(Q, Array("x", "y", "z"))

    def randomPolyQ = ringQ.randomElement(5, 5) + 1

    // some random polynomial composed from some factors
    val polyQ = randomPolyQ * randomPolyQ * randomPolyQ.pow(2)
    // perform square-free factorization
    println(ringQ show FactorSquareFree(polyQ))
    // perform complete factorization
    println(ringQ show Factor(polyQ))


  }

  @Test
  def test24: Unit = {
    import rings.poly.multivar.MultivariateInterpolation._
    import syntax._

    // ring GF(13^6)[x, y, z]
    implicit val ring = MultivariateRing(GF(13, 6, "t"), Array("x", "y", "z"))
    val (x, y, z) = ring("x", "y", "z")

    // coefficient ring GF(13^6)
    val cfRing = ring.cfRing
    val t = cfRing("t")
    // some points for interpolation
    val points: Array[ring.CoefficientType] = {
      // hide implicit ring
      val ring: Any = null
      // enable operations in cfRing
      implicit val _ = cfRing
      Array(1 + t, 2 + t, 3 + t, 12 + t)
    }

    // some values for interpolation
    val values = Array(x + y, x.pow(2) + y * t, y.pow(3), x.pow(4) * t + y)

    // interpolation polynomial values for variable z
    val result = new Interpolation(ring.variable("z"), ring)
      .update(points, values)
      .getInterpolatingPolynomial

    assert(points.zipWithIndex.forall { case (point, i) => result("z" -> point) == values(i) })
  }


  @Test
  def test25: Unit = {
    import syntax._

    // Galois field GF(7^10) represented by univariate polynomials
    // in variable "z" over Z/7 modulo some irreducible polynomial
    // (irreducible polynomial will be generated automatically)
    val gf7_10 = GF(7, 10, "z")
    assert(gf7_10.characteristic == Z(7))
    assert(gf7_10.cardinality == Z(7).pow(10))


    // GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
    val gf7_3 = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
    assert(gf7_3.characteristic == Z(7))
    assert(gf7_3.cardinality == Z(7 * 7 * 7))

    // Mersenne prime 2^107 - 1
    val characteristic = Z(2).pow(107) - 1
    // Galois field GF((2^107 - 1) ^ 16)
    implicit val field = GF(characteristic, 16, "z")

    assert(field.cardinality() == characteristic.pow(16))
  }

  @Test
  def test26: Unit = {
    import syntax._
    {
      // Ring Z/3[x]
      val zp3x = UnivariateRingZp64(3, "x")
      // parse univariate poly from string
      val p1 = zp3x("4 + 8*x + 13*x^2")
      val p2 = zp3x("4 - 8*x + 13*x^2")
      assert(p1 + p2 == zp3x("2 - x^2"))


      // GF(7^3)
      val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
      // GF(7^3)[x]
      val gfx = UnivariateRing(cfRing, "x")
      // parse univariate poly from string
      val r1 = gfx("4 + (8 + z)*x + (13 - z^43)*x^2")
      val r2 = gfx("4 - (8 + z)*x + (13 + z^43)*x^2")
      assert(r1 + r2 == gfx("1 - 2*x^2"))
      val (div, rem) = r1 /% r2
      assert(r1 == r2 * div + rem)
    }
    {
      // Ring Z/3[x, y, z]
      val zp3xyz = MultivariateRingZp64(3, Array("x", "y", "z"))
      // parse univariate poly from string
      val p1 = zp3xyz("4 + 8*x*y + 13*x^2*z^5")
      val p2 = zp3xyz("4 - 8*x*y + 13*x^2*z^5")
      assert(p1 + p2 == zp3xyz("2 - x^2*z^5"))


      // GF(7^3)
      val cfRing = GF(UnivariateRingZp64(7, "t")("1 + 3*t + t^2 + t^3"), "t")
      // GF(7^3)[x, y, z]
      val gfxyz = MultivariateRing(cfRing, Array("x", "y", "z"))
      // parse univariate poly from string
      val r1 = gfxyz("4 + (8 + t)*x*y + (13 - t^43)*x^2*z^5")
      val r2 = gfxyz("4 - (8 + t)*x*y + (13 + t^43)*x^2*z^5")
      assert(r1 + r2 == gfxyz("1 - 2*x^2*z^5"))
      val (div, rem) = r1 /% r2
      assert(r1 == r2 * div + rem)
    }
  }

  @Test
  def test27: Unit = {
    import syntax._

    implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))
    val poly = ring("x^2*y^2 + x + y")

    val rem = poly %% (ring("x + y"), ring("x - y"))

    println(rem)
  }

  @Test
  def test28: Unit = {
    import syntax._

    val inner1 = Zp64(17)
    val inner2 = UnivariateRingZp64(inner1, "t")
    val inner3 = Frac(inner2)
    val inner4 = MultivariateRing(inner3, Array("u, v"))
    implicit val ring = UnivariateRing(inner4, "x")

    val rnd = new Well1024a()
    val poly = (0 to 3).map { _ =>
      val num = inner2.randomElement(33, rnd)
      val den = inner2.randomElement(33, rnd)
      println(inner3 show num)
      println(inner3 show den)
      val rat = new Rational[inner2.ElementType](inner2, num, den)
      println(inner3 show rat)
      val x = ring.randomElement(3, rnd)
      println(ring show x)
      println()
      val r = x.add(inner4.randomElement(3, 3).add(rat))
      r
    }.fold(ring.getOne) { _ * _ }

    println(ring show PolynomialMethods.Factor(poly).size())
  }

  @Test
  def test29: Unit = {
    import syntax._

    // ground field Z/17
    val zpRing = Zp64(17)
    // very large BigInteger
    val n = Z(1).shiftLeft(750)
    // some number
    val r = 19

    // some poly a
    val a = UnivariatePolynomialZp64(1, 2, 3, 4, 3, 2, 1)(zpRing)
    // x^1
    val x = UnivariatePolynomialZp64(0, 1)(zpRing)
    // x - a
    val xa = x - a
    // x^r - 1
    val xr = x.pow(r) - 1

    // Newton iterations for (X^r - 1)^(-1) (used in fast division)
    val xrInv = xr.precomputedInverses

    import UnivariatePolynomialArithmetic._

    // compute (X-a)^n mod (X^r - 1)
    val xan = polyPowMod(a, n, xr, xrInv, true)
    // compute X^n - a mod (X^r - 1)
    val xna = polySubtractMod(polyPowMod(x, n, xr, xrInv, true), a, xr, xrInv, true)

    println(xan)
    println(xna)
  }

  @Test
  def test30: Unit = {
    import syntax._

    /**
      * Some generic function
      *
      * @tparam Poly type of polynomial
      * @tparam Coef type of polynomial coefficients
      */
    def genericFunc[Poly <: IPolynomial[Poly], Coef]
    (poly: Poly)(implicit ring: IPolynomialRing[Poly, Coef]): Poly = {
      // implicit coefficient ring (setups algebraic operators on type Coef)
      implicit val cfRing: Ring[Coef] = ring.cfRing

      // c.c. and l.c. are of type Coef
      val (cc: Coef, lc: Coef) = (poly.cc, poly.lc)
      // do some calculations
      val l: Coef = (cc + lc).pow(10)
      // return the result
      if (ring.characteristic === 2)
        poly - l
      else
        l + l * poly
    }

    implicit val uRing = UnivariateRing(Z, "x")
    // Coef will be inferred as IntZ and Poly as UnivariatePolynomial[IntZ]
    val ur: UnivariatePolynomial[IntZ] = genericFunc(uRing("x^2 + 2"))

    implicit val mRing = MultivariateRingZp64(17, Array("x", "y", "z"))
    // Coef will be inferred as Long and Poly as MultivariatePolynomialZp64
    val mr: MultivariatePolynomialZp64 = genericFunc(mRing("x^2 + 2 + y + z"))
  }


  @Test
  def test31: Unit = {
    import syntax._

    /** Lagrange polynomial interpolation formula */
    def lagrange[Poly <: IUnivariatePolynomial[Poly], Coef]
    (points: Seq[Coef], values: Seq[Coef])
    (implicit ring: IUnivariateRing[Poly, Coef]) = {
      // implicit coefficient ring (setups algebraic operators on type Coef)
      implicit val cfRing: Ring[Coef] = ring.cfRing
      if (!cfRing.isField) throw new IllegalArgumentException
      points.indices
        .foldLeft(ring(0)) { case (sum, i) =>
          sum + points.indices
            .filter(_ != i)
            .foldLeft(ring(values(i))) { case (product, j) =>
              product * (ring.`x` - points(j)) / (points(i) - points(j))
            }
        }
    }

    def interpolate[Poly <: IUnivariatePolynomial[Poly], Coef]
    (points: Seq[(Coef, Coef)])
    (implicit ring: IUnivariateRing[Poly, Coef]) = {
      // implicit coefficient ring (setups algebraic operators on type Coef)
      implicit val cfRing: Ring[Coef] = ring.cfRing
      if (!cfRing.isField) throw new IllegalArgumentException
      points.indices
        .foldLeft(ring(0)) { case (sum, i) =>
          sum + points.indices
            .filter(_ != i)
            .foldLeft(ring(points(i)._2)) { case (product, j) =>
              product * (ring.`x` - points(j)._1) / (points(i)._1 - points(j)._1)
            }
        }
    }

    // coefficient ring Frac(Z/13[a,b,c])
    val cfRing = Frac(MultivariateRingZp64(2, Array("a", "b", "c")))
    val (a, b, c) = cfRing("a", "b", "c")

    implicit val ring = UnivariateRing(cfRing, "x")
    // interpolate with Lagrange formula
    val data = Seq(a -> b, b -> c, c -> a)
    val poly = interpolate(data)
    assert(data.forall { case (p, v) => poly.eval(p) == v })

    println(poly)


    //UnivariatePolynomial[Rational[MultivariatePolynomialZp64]]

    //    // Q[x]
    //    implicit val ringA = UnivariateRing(Q, "x")
    //    // Coef will be inferred as IntZ and Poly as UnivariatePolynomial[IntZ]
    //    val ur: UnivariatePolynomial[Rational[IntZ]]
    //    = lagrange(
    //      Array(1, 2, 3).map(Q(_)),
    //      Array(5, 4, 3).map(Q(_))
    //    )
  }

  @Test
  def test32: Unit = {
    import syntax._
    implicit val ring = MultivariateRing(Zp(2), Array("x", "y", "z"))
    val poly = ring("x^2*y - z^2*x^3 - y^2*z - 1").pow(20) - 1
    println(PolynomialMethods.Factor(poly))
  }

  @Test
  def test33: Unit = {
    import syntax._

    implicit val ring = Frac(MultivariateRing(Zp(2), Array("x", "y", "z")))
    val (x, y, z) = ring("x", "y", "z")
    val lhs =
      Array(
        Array(x + z, y * z, z - x * y),
        Array(x - y, x / y, z + x / y),
        Array(x * y, x + y, z / x + y)
      )
    val rhs = Array(x, y, z)

    val solution = LinearSolver.solve[ring.ElementType](ring, lhs, rhs)

    println(solution.mkString("\n"))
  }

  @Test
  def test34: Unit = {
    import syntax._

    implicit val ring = UnivariateRing(Frac(MultivariateRing(Zp(47), Array("x", "y", "z"))), "t")
    val (x, y, z, t) = ring("(x + y)", "(y)", "(z)", "t")
    println(ring show IrreduciblePolynomials.randomIrreduciblePolynomial[ring.cfRing.ElementType](ring.cfRing, 4, new Random()))
    val s = rings.ChineseRemainders.ChineseRemainders[ring.ElementType](ring, x + y + t.pow(2), y + t, x + t, z + t.pow(2))
    println(s)
  }

  @Test
  def test35: Unit = {
    import syntax._

    val cfRing = Zp(17)
    implicit val ring = UnivariateRing(cfRing, "x")


    val random = new Random
    val nextPrime = (degree: Int) => {
      import IrreduciblePolynomials._
      randomIrreduciblePolynomial[cfRing.ElementType](cfRing, degree, random)
    }

    val primes = Array(nextPrime(5), nextPrime(6), nextPrime(7))
    val remainders = primes.map(_.degree()).map(i => ring.randomElement(random.nextInt(i), random))

    import rings.ChineseRemainders._

    val crt = ChineseRemainders[ring.ElementType](ring, primes, remainders)
    assert(crt.degree() == primes.map(_.degree()).sum - 1)
    assert((primes zip remainders).forall { case (prime, rem) => crt % prime == rem })
  }

  @Test
  def test36: Unit = {
    import syntax._

    // Galois field %$\color{commcolor}\GF(17, 3)$%
    implicit val gf = GF(17, 3, "t")
    // parse ring elements from their string representation
    val (t, t1): (UnivariatePolynomialZp64, UnivariatePolynomialZp64) = gf("t", "1 + t^2")
    // do some basic math (+-*/)
    val gfOps = 1 + t1 - t.pow(22) / (1 + t + t1.pow(999))

    // multivariate ring %$\color{commcolor}\GF(17, 3)[x, y, z]$%
    implicit val ring = MultivariateRing(gf, Array("x", "y", "z"), LEX)
    // parse some polynomial from string
    val p1 = ring("(1 + t + t^2)*x*y^2*z + (1 + t)*x^5*z*y^6 + 1")
    val (x, y, z) = ring("x", "y", "z")

    // do some basic math (+-*/) with elements of ring
    val p2 = p1.pow(2) + (t + 1) * x.pow(2) * y.pow(2) + (t.pow(9) + 1) * z.pow(7)


    val (div, rem) = p2 /% p1

    val ideal = (x + y + z, x - y - z, y.pow(2) - z.pow(2))

    val (div1, div2, div3, rem1) = p2 /%/% ideal
    assert(div1 * ideal._1 + div2 * ideal._2 + div3 * ideal._3 + rem1 == p2)

    val p3 = p2 %% ideal
    assert(p3 == rem1)

    val gcd = PolynomialMethods.PolynomialGCD(p1 * p3, p2 * p3)
    assert(gcd % p3 === 0)

    val gcd2 = PolynomialMethods.PolynomialGCD(p1 * p3, p2 * p3 + 1)
    assert(gcd2.isConstant)

    val poly = p1 * p2.pow(2) * p3.pow(3)
    println(poly.size())

    println(poly.degrees().mkString(","))
    //    val factors = PolynomialMethods.Factor(poly)


    implicit val ratFuncs = Frac(ring)
    // x, y, z and t as elements of Frac[GF(17, 3)[x,y,z]]
    val (rx, ry, rz, rt) = ratFuncs(x, y, z, ring(t))


    val lhs =
      Array(
        Array(rt + rx + rz, ry * rz, rz - rx * ry),
        Array(rx - ry - rt, rx / ry, rz + rx / ry),
        Array(rx * ry / rt, rx + ry, rz / rx + ry)
      )
    val rhs = Array(rx, ry, rz)
    val solution = LinearSolver.solve[ratFuncs.ElementType](ratFuncs, lhs, rhs)


    def solveDiophantine[E](fi: Seq[E])(implicit ring: Ring[E]) =
      fi.foldLeft((ring(0), Seq.empty[E])) { case ((gcd, seq), f) =>
        val xgcd = ring.extendedGCD(gcd, f)
        (xgcd(0), seq.map(_ * xgcd(1)) :+ xgcd(2))
      }

    def apart[E](frac: Rational[E]) = {
      implicit val ring: Ring[E] = frac.ring
      val factors = ring.factor(frac.denominator).map { case (f, exp) => f.pow(exp) }
      val (gcd, nums) = solveDiophantine(factors.map(frac.denominator / _))
      val (ints, rats) = (nums zip factors)
        .map { case (num, den) => Rational(frac.numerator * num, den * gcd) }
        .flatMap(_.normal) // reduce Rational to normal form ( /% )
        .partition(_.isIntegral) // select fractions and integers
      rats :+ ints.foldLeft(Rational(ring(0)))(_ + _)
    }


    implicit val uRing = UnivariateRing(ratFuncs, "W")
    assert(uRing.isEuclideanRing)
    val W = uRing("W")


    val wpoly = W.pow(2) + 1
    val fracs = apart(Rational(W + 1, wpoly))
    println(fracs)

    // univariate polynomial rx/ry
    val prime1 = rx / ry + (rt / rx) * W + W.pow(2) + (rz + 1) * W.pow(3)
    val prime2 = prime1 * W.pow(3) + 1

    println(apart(Rational(W, prime1 * prime2)))


    import IrreduciblePolynomials._
    assert(irreducibleQ(prime1) && irreducibleQ(prime2))

    val urem1 = rx + ry + rt * W.pow(2)
    val urem2 = rx + rz + W

    val primes = Array(prime1, prime2)
    val rems = Array(urem1, urem2)

    import rings.ChineseRemainders._

    val crt = ChineseRemainders(uRing, primes, rems)
    println(uRing show crt)

    assert(crt.degree() == primes.map(_.degree()).sum - 1)
    assert((primes zip rems).forall { case (prime, re) => crt % prime == re })

    import MultivariateInterpolation._

    val interpolation = new Interpolation(0, ring)
    interpolation.update(t, y)
    interpolation.update(t + 1, y + z)
    interpolation.update(t + 2, y.pow(2) + z.pow(2))
    val result = interpolation.getInterpolatingPolynomial
    assert(result("x" -> t) == y)
    assert(result("x" -> (t + 1)) == y + z)

    ring.isEuclideanRing // whether ring is Euclidean ring
    ring.isField // whether ring is a field
    ring.isFinite // whether ring is finite
    ring.cardinality // ring cardinality
    ring.characteristic // ring characteristic

  }


  @Test
  def test37: Unit = {
    import syntax._
    import PolynomialMethods._
    import scala.collection.JavaConverters._

    def solveDiophantine[E](fi: Seq[E])(implicit ring: Ring[E]) =
      fi.foldLeft((ring(0), Seq.empty[E])) { case ((gcd, seq), f) =>
        val xgcd = ring.extendedGCD(gcd, f)
        (xgcd(0), seq.map(_ * xgcd(1)) :+ xgcd(2))
      }


    def apart[E](frac: Rational[E]) = {
      implicit val ring: Ring[E] = frac.ring
      val factors = ring.factor(frac.denominator).map { case (f, exp) => f.pow(exp) }
      val (gcd, nums) = solveDiophantine(factors.map(frac.denominator / _))
      val (ints, rats) = (nums zip factors)
        .map { case (num, den) => Rational(frac.numerator * num, den * gcd) }
        .flatMap(_.normal) // reduce Rational to normal form ( /% )
        .partition(_.isIntegral) // select fractions and integers
      rats :+ ints.foldLeft(Rational(ring(0)))(_ + _)
    }

    val qFracs = apart(Q("1234213/2341352"))

    val ring = Frac(UnivariateRingZp64(17, "x"))
    val pFracs = apart(ring("1/(3 - 3*x^2 - x^3 + x^5)"))


    println(qFracs)
    println(pFracs)
  }
}