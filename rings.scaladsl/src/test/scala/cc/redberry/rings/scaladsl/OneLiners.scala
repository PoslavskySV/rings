package cc.redberry.rings.scaladsl

import cc.redberry.rings.Rational
import cc.redberry.rings.poly.PolynomialMethods.PolynomialGCD
import cc.redberry.rings.poly.Util
import org.junit.Test

/**
  *
  * @author Stanislav Poslavsky
  * @since 1.0
  */
class OneLiners {

  import syntax._


  @Test
  def test1: Unit = {
    def reconstruct[Poly <: IUnivariatePolynomial[Poly]]
    (n: Poly, modulus: Poly, numBound: Int, denBound: Int)(implicit ring: Ring[Poly]): Option[Rational[Poly]] = {
      var v = (modulus, ring(0))
      var w = (n, ring(1))
      while (w._1.degree > numBound) {
        val q = (v._1 /% w._1)._1
        val z = (v._1 - (q * w._1), v._2 - (q * w._2))
        v = w
        w = z
      }
      if (w._2.degree <= denBound && PolynomialGCD(w._1, w._2).isConstant)
        Some(Rational(w._1, w._2))
      else
        None
    }

    def padeApproximant[Poly <: IUnivariatePolynomial[Poly], Cf](poly: Poly, point: Cf, numBound: Int, denBound: Int)
                                                                (implicit ring: IUnivariateRing[Poly, Cf]): Option[Rational[Poly]] = {
      val mod: Poly = (ring.`x` - point).pow(numBound + denBound + 1)
      reconstruct(poly, mod, numBound, denBound)
    }


    implicit val ring = UnivariateRingZp64(17, "x")

    val poly = ring("1 + x + 2*x^2 + x^3 + x^4")

    val pade = padeApproximant(poly, ring.cfRing(0), 3, 5)
    println(pade)
  }

  /**
    * Solves equation %$\color{commcolor} \sum f_i s_i  = gcd(f_1, \dots, f_N)$% for given %\color{commcolor}$f_i$% and unknown %\color{commcolor}$s_i$%
    *
    * @return a tuple (gcd, solution)
    */
  def solveDiophantine[E](fi: Seq[E])(implicit ring: Ring[E]) =
    fi.foldLeft((ring(0), Seq.empty[E])) { case ((gcd, seq), f) =>
      val xgcd = ring.extendedGCD(gcd, f)
      (xgcd(0), seq.map(_ * xgcd(1)) :+ xgcd(2))
    }

  /** Computes partial fraction decomposition of given rational */
  def apart[E](frac: Rational[E]) = {
    implicit val ring: Ring[E] = frac.ring
    val factors = ring.factor(frac.denominator).map { case (f, exp) => f.pow(exp) }
    val (gcd, nums) = solveDiophantine(factors.map(frac.denominator / _))
    val (ints, rats) = (nums zip factors)
      .map { case (num, den) => Rational(frac.numerator * num, den * gcd) }
      .flatMap(_.normal) // extract integral parts from fractions
      .partition(_.isIntegral) // separate integrals and fractions
    rats :+ ints.foldLeft(Rational(ring(0)))(_ + _)
  }

  def apart[Term <: AMonomial[Term], Poly <: AMultivariatePolynomial[Term, Poly], Cf]
  (frac: Rational[Poly], variable: Int)(implicit ring: IMultivariateRing[Term, Poly, Cf]): Seq[Rational[Poly]] = {

    import Conversions._

    val num = frac.numerator
    val den = frac.denominator

    implicit val uRingZ: UnivariateRing[Poly] = asUnivariate[Poly, Cf](ring, variable)
    implicit val cfField: Ring[Rational[Poly]] = Frac(uRingZ.cfRing)
    implicit val uRingQ: UnivariateRing[Rational[Poly]] = uRingZ.copy(cfRing = cfField)

    val uNum: UnivariatePolynomial[Poly] = asUnivariate(num, variable)
    val uDen: UnivariatePolynomial[Poly] = asUnivariate(den, variable)

    val uNumF: UnivariatePolynomial[Rational[Poly]] = uNum.map(cfField, p => Rational(p)(uRingZ.cfRing))
    val uDenF: UnivariatePolynomial[Rational[Poly]] = uDen.map(cfField, p => Rational(p)(uRingZ.cfRing))

    val seq = apart(Rational(uNumF, uDenF))
    seq
      .map { frac =>
        val fnum: UnivariatePolynomial[Rational[Poly]] = frac.numerator
        val fden: UnivariatePolynomial[Rational[Poly]] = frac.denominator

        val fnumc = Util.toCommonDenominator(fnum)
        val fdenc = Util.toCommonDenominator(fden)

        val rnum = fromUnivariate(fnumc._1 * fdenc._2, variable)
        val rden = fromUnivariate(fdenc._1 * fnumc._2, variable)

        Rational(rnum, rden)
      }
  }

  @Test
  def test3: Unit = {
    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))
    val (x, y, z) = ring("x", "y", "z")

    val den = (x.pow(22) * y - z - y) * (x.pow(15) + y.pow(2) * z + 1)
    val initial = Rational(1, den)
    for (i <- 1 to 10)
      println(util.timing({ apart[ring.MonomialType, ring.ElementType, ring.CoefficientType](initial, 0) })._1)
    val (t, fracs) = util.timing { apart[ring.MonomialType, ring.ElementType, ring.CoefficientType](initial, 0) }

    println(t)

    println(fracs)

    implicit val field = Frac(ring)
    println(initial == fracs.foldLeft(field(0))(_ + _))
  }

  @Test
  def test2: Unit = {
    /**
      * Solves equation %$\color{commcolor} \sum f_i s_i  = gcd(f_1, \dots, f_N)$% for given %\color{commcolor}$f_i$% and unknown %\color{commcolor}$s_i$%
      *
      * @return a tuple (gcd, solution)
      */
    def solveDiophantine[E](fi: Seq[E])(implicit ring: Ring[E]) =
      fi.foldLeft((ring(0), Seq.empty[E])) { case ((gcd, seq), f) =>
        val xgcd = ring.extendedGCD(gcd, f)
        (xgcd(0), seq.map(_ * xgcd(1)) :+ xgcd(2))
      }

    /** Computes partial fraction decomposition of given rational */
    def apart[E](frac: Rational[E]) = {
      implicit val ring: Ring[E] = frac.ring
      val factors = ring.factor(frac.denominator).map { case (f, exp) => f.pow(exp) }
      val (gcd, nums) = solveDiophantine(factors.map(frac.denominator / _))
      val (ints, rats) = (nums zip factors)
        .map { case (num, den) => Rational(frac.numerator * num, den * gcd) }
        .flatMap(_.normal) // extract integral parts from fractions
        .partition(_.isIntegral) // separate integrals and fractions
      rats :+ ints.foldLeft(Rational(ring(0)))(_ + _)
    }


    //    // partial fraction decomposition for rationals
    //    // gives List(184/479, (-10)/13, 1/8, (-10)/47, 1)
    //    val qFracs = apart(Q("1234213 / 2341352"))
    //
    //    // partial fraction decomposition for rational functions
    //    val ufRing = Frac(UnivariateRingZp64(17, "x"))
    //    // gives List(4/(16+x), 1/(10+x), 15/(1+x), (14*x)/(15+7*x+x^2))
    //    val pFracs = apart(ufRing("1 / (3 - 3*x^2 - x^3 + x^5)"))


    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))
    val (x, y, z) = ring("x", "y", "z")

    val den = (x.pow(2) * y - z - y) * (x + y.pow(2) * z + 1)
    import Conversions._

    implicit val uRing: UnivariateRing[MultivariatePolynomial[Rational[IntZ]]] = asUnivariate(ring, 0)
    implicit val uFracs: UnivariateRing[Rational[MultivariatePolynomial[Rational[IntZ]]]] = uRing.copy(cfRing = Frac(uRing.cfRing))
    val pp: UnivariatePolynomial[MultivariatePolynomial[Rational[IntZ]]] = asUnivariate(den, 0)
    val fff: Ring[Rational[MultivariatePolynomial[Rational[IntZ]]]] = uFracs.cfRing
    val uPoly = pp.map(fff, f => Rational(f)(uRing.cfRing))


    println(uFracs show uPoly)
    val seq = apart(Rational(1, uPoly))
    println(seq)
    //    implicit val uRing = asUnivariate(ring, 0)
    //
    //    val uPoly = asUnivariate(den, 0).mapCoefficients(uRing.co)
    //
    //    val rat = Rational(1, uPoly)(uRing)
    //    println(apart(rat))

  }
}
