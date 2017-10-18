.. image:: https://circleci.com/gh/PoslavskySV/rings.svg?style=svg
   :target: https://circleci.com/gh/PoslavskySV/rings

.. role:: feature
  :class: feature

.. |br| raw:: html

   <br/>

.. |Rings| raw:: html

   <span style="color:#003366;font: normal Lucida Sans, Lucida Sans Unicode, Lucida Grande, sans-serif; font-variant:small-caps;" >Rings</span>

.. |_| unicode:: 0xA0 
   :trim:

.. |____| replace:: |_|


Rings: efficient Java/Scala library for polynomial rings
########################################################

|Rings| is an efficient implementation of univariate and multivariate polynomial algebra over arbitrary coefficient rings. It make use of asymptotically fast algorithms for basic algebraic operations as well as for advanced methods like GCDs and polynomial factorization. Performance achieved in |Rings| is comparable to well known solutions like Singular/NTL/FLINT/Maple/Mathematica.

The key features of |Rings| include:

 * `Polynomials → <http://rings.readthedocs.io/en/latest/guide.html#ref-basics-polynomials>`_ |br| :feature:`Univariate and multivariate polynomials over arbitrary coefficient rings`

 * `Polynomial GCD → <http://rings.readthedocs.io/en/latest/guide.html#ref-polynomial-methods>`_ |br| :feature:`Polynomial GCD over arbitrary coefficient domains`

 * `Univariate factorization → <http://rings.readthedocs.io/en/latest/guide.html#ref-univariate-factorization>`_ |br| :feature:`Univariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
   
 * `Multivariate factorization → <http://rings.readthedocs.io/en/latest/guide.html#ref-multivariate-factorization>`_ |br| :feature:`Multivariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
 
 * `Algebra → <http://rings.readthedocs.io/en/latest/guide.html#ref-rings>`_ |br| :feature:`Arbitrary rings, Galois fields etc`

 * `Scala DSL → <http://rings.readthedocs.io/en/latest/guide.html#ref-scala-dsl>`_ |br| :feature:`Powerful domain specific language in Scala`
   
 * `Fast → <http://rings.readthedocs.io/en/latest/quickstart.html#ref-some-benchamrks>`_ |br| :feature:`Really fast library suitable for real-world computational challenges`


The full documentation is available at `http://rings.readthedocs.io <https://rings.readthedocs.io>`_. 


======
Set up
======


Interactive Rings shell
=======================

To taste what |Rings| can do, one can try interactive |Rings| session with `Ammonite REPL <http://ammonite.io>`_. Type the following commands at the prompt to install |Rings|\ *.repl*:

.. code-block:: bash

    $ sudo curl -L -o /usr/local/bin/amm https://git.io/v5Tct && sudo chmod +x /usr/local/bin/amm
    $ sudo curl -L -o /usr/local/bin/rings.repl https://git.io/vdQ6P && chmod +x /usr/local/bin/rings.repl

and run:

.. code-block:: scala

    $ rings.repl
    Loading...
    Rings 2.0: efficient Java/Scala library for polynomial rings

    @ implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))
    ring: MultivariateRing[IntZ] = MultivariateRing(Z, Array("x", "y", "z"), LEX)

    @ val poly1 = ring("x + y - z").pow(8) 
    poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*y^4*z^4-56...

    @ val poly2 = ring("x - y + z").pow(8) 
    poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*y^4*z^4-56...

    @ Factor(poly1 - poly2)
    res13: FactorDecomposition[MultivariatePolynomial[IntZ]] = 16*x*(-z+y)*(z^2-2*y*z+y^2+x^2)*(z^4-4*y*z^3+6*y^2*z^2-4*y^3*z+y^4+6*x^2*z^2-12*x^2*y*z+6*x^2*y^2+x^4)



Java/Scala library
==================

|Rings| are currently available for Java and Scala. To get started with Scala SBT, simply add the following dependence to your ``build.sbt`` file:

.. code-block:: scala

    libraryDependencies += "cc.redberry" % "rings.scaladsl" % "2.0"

For using |Rings| solely in Java there is Maven artifact:

.. code-block:: xml

    <dependency>
        <groupId>cc.redberry</groupId>
        <artifactId>rings</artifactId>
        <version>2.0</version>
    </dependency>



======================================
Examples: algebra, GCDs, factorization
======================================

Below examples can be evaluated directly in the |Rings|\ *.repl*. If using in Scala, the following preambula will import all required things from |Rings| library:

.. code-block:: scala

    import cc.redberry.rings

    import rings.poly.PolynomialMethods._
    import rings.scaladsl._
    import syntax._


Java examples can be found in the `complete documentation pages <https://rings.readthedocs.io>`_. 

----

Do some algebra in Galois field :math:`GF(17^{9})`:

.. code-block:: scala

    // GF(17^9) (irreducible poly in Z/17[x] will be generated automaticaly)
    implicit val ring = GF(17, 9, "x")

    // some random element from ring
    val a = ring.randomElement()
    val b = a.pow(1000)
    val c = 1 / b

    assert ( b * c === 1)

    // explicitly parse ring element from string
    val d = ring("1 + x + x^2 + x^3 + 15*x^999")
    // do some math ops
    val some = a / (b + c) + a.pow(6) - a * b * c * d


----

Some math with multivariate polynomials from :math:`Z[x, y, z]`:

.. code-block:: scala

    // Z[x, y, z]
    implicit val ring = MultivariateRing(Z, Array("x", "y", "z")) 

    val (x, y, z) = ring("x", "y", "z") 

    // do some math
    val a = (x + y + z).pow(2) - 1 
    val b = (x - y - z - 1).pow(2) + x + y + z - 1 
    val c = (a + b + 1).pow(9) - a - b - 1

    // reduce c modulo a and b (multivariate division with remainder)
    val (div1, div2, rem) = c /%/% (a, b)


----

Univariate extended GCD in :math:`Z_{17}[x]`:

.. code-block:: scala

    // ring Z/17[x]
    implicit val ring = UnivariateRingZp64(17, "x")

    val x = ring("x")
    
    val poly1 = 1 + x + x.pow(2) + x.pow(3)
    val poly2 = 1 + 2 * x + 9 * x.pow(2)
    val (gcd, s, t) = PolynomialExtendedGCD(poly1, poly2).tuple3

    println((gcd, s, t))

----

Multivariate GCD in :math:`Z[a, b, c]`:

.. code-block:: scala

    // ring Z[a, b, c]
    implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))

    val poly1 = ring("-b-b*c-b^2+a+a*c+a^2")
    val poly2 = ring("b^2+b^2*c+b^3+a*b^2+a^2+a^2*c+a^2*b+a^3")

    val gcd   = PolynomialGCD(poly1, poly2)

    println(s"gcd: ${ring show gcd}")


----

Factor polynomial in :math:`Z_{17}[x]`:

.. code-block:: scala

    // ring Z/17[x]
    implicit val ring = UnivariateRingZp64(17, "x")x

    val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 - x^6 + 10*x^7 + x^8")

    // factorize poly
    val factors = Factor(poly)

    println(factors)


Coefficient rings with arbitrary large characteristic are available:


.. code-block:: scala

    // coefficient ring Z/1237940039285380274899124357 (the next prime to 2^100)
    val modulus = Z("1267650600228229401496703205653")
    val cfRing  = Zp(modulus)

    // ring Z/1237940039285380274899124357[x]
    implicit val ring = UnivariateRing(cfRing, "x")

    val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 + 16*x^6 + 27*x^7 + 18*x^8")
    
    println(Factor(poly))



(large primes can be generated with ``BigPrimes.nextPrime`` method, see `Prime numbers <http://rings.readthedocs.io/en/latest/guide.html#ref-primes>`_).


----

Factor polynomial in :math:`Z_{2}[x, y, z]`:

.. code-block:: scala

    // ring Z/2[x, y, z]
    implicit val ring = MultivariateRingZp64(2, Array("x", "y", "z"))

    val (x, y, z) = ring("x", "y", "z")
    
    val factors = Factor(1 + (1 + x + y + z).pow(2) + (x + y + z).pow(4))

    println(factors)

----

Factor polynomial in :math:`Z[a, b, c]`:

.. code-block:: scala

    // ring Z[a, b, c]
    implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))

    val (a, b, c) = ring("a", "b", "c")
    
    val factors = Factor(1 - (1 + a + b + c).pow(2) - (2 + a + b + c).pow(3))

    println(ring show factors)


----

Factor polynomial in :math:`Q[x, y, z]`:

.. code-block:: scala

    // ring Q[x, y, z]
    implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))

    val poly = ring(
      """
        |(1/6)*y*z + (1/6)*y^3*z^2 - (1/2)*y^6*z^5 - (1/2)*y^8*z^6
        |-(1/3)*x*z - (1/3)*x*y^2*z^2 + x*y^5*z^5 + x*y^7*z^6
        |+(1/9)*x^2*y^2*z - (1/3)*x^2*y^7*z^5 - (2/9)*x^3*y*z
        |+(2/3)*x^3*y^6*z^5 - (1/2)*x^6*y - (1/2)*x^6*y^3*z
        |+x^7 + x^7*y^2*z - (1/3)*x^8*y^2 + (2/3)*x^9*y
      """.stripMargin)

    val factors = Factor(poly)

    println(factors)


----

Polynomial rings over :math:`Z` and :math:`Q`:

.. code-block:: scala

    // Ring Z[x]
    UnivariateRing(Z, "x")
    // Ring Z[x, y, z]
    MultivariateRing(Z, Array("x", "y", "z"))
    // Ring Q[a, b, c]
    MultivariateRing(Q, Array("a", "b", "c"))

Polynomial rings over :math:`Z_p`:

.. code-block:: scala

    // Ring Z/3[x]
    UnivariateRingZp64(3, "x")
    // Ring Z/3[x, y, z]
    MultivariateRingZp64(3, Array("x", "y", "z"))
    // Ring Z/p[x, y, z] with p = 2^107 - 1 (Mersenne prime)
    MultivariateRing(Zp(Z(2).pow(107) - 1), Array("x", "y", "z"))


Galois fields:

.. code-block:: scala

    // Galois field with cardinality 7^10 
    // (irreducible polynomial will be generated automatically)
    GF(7, 10, "x")
    // GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
    GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")


Fractional fields:

.. code-block:: scala

    // Field of fractions of univariate polynomials Z[x]
    Frac(UnivariateRing(Z, "x"))
    // Field of fractions of multivariate polynomials Z/19[x, y, z]
    Frac(MultivariateRingZp64(19, Array("x", "y", "z")))

----

Ring of univariate polynomials over elements of Galois field :math:`GF(7^{3})[x]`:

.. code-block:: scala

    // Elements of GF(7^3) are represented as polynomials
    // over "z" modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
    val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")

    assert(cfRing.characteristic().intValue() == 7)
    assert(cfRing.cardinality().intValue() == 343)

    // Ring GF(7^3)[x]
    implicit val ring = UnivariateRing(cfRing, "x")

    // Coefficients of polynomials in GF(7^3)[x] are elements of GF(7^3)
    val poly = ring("1 - (1 - z^3) * x^6 + (1 - 2*z) * x^33 + x^66")

    // factorize poly (in this examples there will be 9 factors)
    val factors = Factor(poly)
    println(s"${ring show factors}")

----

Ring of multivariate polynomials over elements of Galois field :math:`GF(7^{3})[x, y, z]`:

.. code-block:: scala

    // Elements of GF(7^3) are represented as polynomials
    // over "z" modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
    val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
    // Ring GF(7^3)[x]
    implicit val ring = MultivariateRing(cfRing, Array("a", "b", "c"))

    // Coefficients of polynomials in GF(7^3)[x] are elements of GF(7^3)
    val poly = ring("1 - (1 - z^3) * a^6*b + (1 - 2*z) * c^33 + a^66")

.. _ref-some-benchamrks:

===============
Some benchmarks
===============

In the following plots performance of |Rings| is compared to Wolfram Mathematica 11. All tests were performed on MacBook Pro (15-inch, 2017), 3,1 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3. The code of benchmarks can be found on GitHub. In all benchamrks random polynomials were used.


.. figure:: doc/_static/bench_gcd_Z.png
   :scale: 50%
   :align: center

   Polynomial GCD performance on random sparse multivariate polynomials in :math:`Z[x, y, z]` with about 100 terms, and degree equal to 20 in each variable. |Rings| is about 3 times faster.


.. figure:: doc/_static/bench_gcd_Z2.png
   :scale: 50%
   :align: center

   Polynomial GCD performance on random sparse multivariate polynomials in :math:`Z_2[x, y, z]` with about 100 terms, and degree equal to 20 in each variable. Points marked with red color are those where Mathematica failed to obtain result in less than 3 minutes. |Rings| is about 100 times faster.


.. figure:: doc/_static/bench_fac_uni_Zp.png
   :scale: 50%
   :align: center

   Univariate factorization performance on random polynomials in :math:`Z_{32771}[x]` of degree 250. |Rings| is about 15% slower (this difference remains the same for polynomials of larger degrees).


.. figure:: doc/_static/bench_fac_multi_Z.png
   :scale: 50%
   :align: center

   Multivariate factorization performance on random sparse polynomials in :math:`Z[x_1, x_2, x_3, x_4]` with at least 2 factors with size 100 and degree 10 in each variable.  |Rings| is about 9 times faster.


========================================
Index of algorithms implemented in Rings
========================================



Univariate polynomials
======================

1. *Karatsuba multiplication* |____| (Sec. 8.1 in [GaGe03]_) |br| used with some adaptations for multiplication of univariate polynomials: 

 - `UnivariatePolynomial.multiply <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariatePolynomial.java>`_
 - `UnivariatePolynomialZp64.multiply <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariatePolynomialZp64.java>`_
     
2. *Half-GCD and Extended Half-GCD* |____| (Sec. 11 in [GaGe03]_) |br| used with adaptations inspired by [ShoNTL]_ implementation for univariate GCD:

 - `UnivariateGCD.HalfGCD  <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_
 - `UnivariateGCD.ExtendedHalfGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_
 
3. *Subresultant polynomial remainder sequences* |____| (Sec. 7.3 in [GeCL92]_):

 - `UnivariateGCD.SubresultantRemainders <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_

4. *Modular GCD in* :math:`Z[x]` *and* :math:`Q[x]` |____| (Sec. 6.7 in [GaGe03]_, small primes version):

 - `UnivariateGCD.ModularGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_

5. *Fast univariate division with Newton iteration* |____| (Sec. 9.1 in [GaGe03]_) |br| used everywhere where multiple divisions (remainders) by the same divider are performed:

 - `UnivariateDivision.fastDivisionPreConditioning <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateDivision.java>`_
 - `UnivariateDivision.divideAndRemainderFast <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateDivision.java>`_
 
6. *Univariate square-free factorization in zero characteristic (Yun's algorithm)* |____| (Sec. 14.6 in [GaGe03]_):

 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
     
7. *Univariate square-free factorization in non-zero characteristic (Musser's algorithm)* |____| (Sec. 8.3 in [GeCL92]_, [Muss71]_):

 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationMusser <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
 
8. *Distinct-degree factorization* |____| (Sec. 14.2 in [GaGe03]_) |br| plain version and adapted version with precomputed :math:`x`-powers (used by default):

 - `DistinctDegreeFactorization.DistinctDegreeFactorizationPlain <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_
 - `DistinctDegreeFactorization.DistinctDegreeFactorizationPrecomputedExponents <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_

9. *Shoup's baby-step giant-step algorithm for distinct-degree factorization* |____| ([Shou95]_) |br| used for factorization over fields with large cardinality:

 - `DistinctDegreeFactorization.DistinctDegreeFactorizationShoup <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_

10. *Univariate modular composition* |br| plain algorithm with Horner schema:
 
 - `ModularComposition.compositionHorner <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/ModularComposition.java>`_

11. *Brent-Kung univariate modular composition* |____| ([BreK98]_, [Shou95]_):

 - `ModularComposition.compositionBrentKung <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/ModularComposition.java>`_

12. *Cantor-Zassenhaus algorithm (equal-degree splitting)* |____| (Sec. 14.3 in [GaGe03]_) |br| both for odd and even characteristic:

 - `EqualDegreeFactorization.CantorZassenhaus <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/EqualDegreeFactorization.java>`_

13. *Univaraite linear p-adic Hensel lifting* |____| (Sec. 6.5 in [GeCL92]_):

 - `univar.HenselLifting.createLinearLift <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_
 - `univar.HenselLifting.liftFactorization <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_

14. *Univaraite quadratic p-adic Hensel lifting* |____| (Sec. 15.4-15.5 in [GaGe03]_):

 - `univar.HenselLifting.createQuadraticLift <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_
 - `univar.HenselLifting.liftFactorization <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_

15. *Univariate polynomial factorization over finite fields* |br| uses Musser's square free factorization followed by distinct-degree factorization (either :math:`x`-powers or Shoup's algorithm) followed by Cantor-Zassenhaus equal-degree factorization:

 - `UnivariateFactorization.FactorInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_

16. *Univariate polynomial factorization over Z and Q* |br| uses factorization modulo small prime followed by Hensel lifting (adaptive linear/quadratic) and naive recombination:

 - `UnivariateFactorization.FactorInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_
 - `UnivariateFactorization.FactorInQ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_

17. *Univariate irreducibility test* |____| (Sec. 14.9 in [GaGe03]_):

 - `IrreduciblePolynomials.irreducibleQ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/IrreduciblePolynomials.java>`_

18. *Ben-Or’s generation of irreducible polynomials* |____| (Sec. 14.9 in [GaGe03]_):

 - `IrreduciblePolynomials.randomIrreduciblePolynomial <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/IrreduciblePolynomials.java>`_

19. *Univariate polynomial interpolation* |br| Lagrange and Newton methods:

 - `UnivariateInterpolation <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateInterpolation.java>`_


Multivariate polynomials
========================


20. *Brown GCD over finite fields* |____| ([Brow71]_, Sec. 7.4 in [GeCL92]_, [Yang09]_):

 - `MultivariateGCD.BrownGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

21. *Zippel's sparse GCD over finite fields* |____| ([Zipp79]_, [Zipp93]_, [dKMW05]_, [Yang09]_) |br| both for monic (with fast Vandermonde systems) and non-monic (LINZIP) cases:

 - `MultivariateGCD.ZippelGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

22. *Extended Zassenhaus GCD (EZ-GCD) over finite fields* |____| (Sec. 7.6 in [GeCL92]_, [MosY73]_):

 - `MultivariateGCD.EZGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

23. *Enhanced Extended Zassenhaus GCD (EEZ-GCD) over finite fields* |____| ([Wang80]_):

 - `MultivariateGCD.EEZGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

24. *Modular GCD over Z with sparse interpolation* |____| ([Zipp79]_, [Zipp93]_, [dKMW05]_) |br| (the same interpolation techniques as in `ZippelGCD` is used):

 - `MultivariateGCD.ModularGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

25. *Kaltofen's & Monagan's generic modular GCD* |____| ([KalM99]_) |br| used for computing multivariate GCD over finite fields of very small cardinality:

 - `MultivariateGCD.ModularGCDInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

26. *Multivariate square-free factorization in zero characteristic (Yun's algorithm)* |____| ([LeeM13]_):

 - `MultivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_

27. *Multivariate square-free factorization in non-zero characteristic (Musser's algorithm)* |____| ([Muss71]_, Sec. 8.3 in [GeCL92]_):

- `MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_
- `MultivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_

28. *Bernardin's fast dense multivariate Hensel lifting* |____| ([Bern99]_, [LeeM13]_) |br| both for bivariate case (original Bernardin's paper) and multivariate case (Lee thesis) and both with and without precomputed leading coefficients:

- `multivar.HenselLifting <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/HenselLifting.java>`_

29. *Fast dense bivariate factorization with recombination* |____| ([Bern99]_, [LeeM13]_):

- `MultivariateFactorization.bivariateDenseFactorSquareFreeInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_
- `MultivariateFactorization.bivariateDenseFactorSquareFreeInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

30. *Kaltofen's multivariate factorization in finite fields* |____| ([Kalt85]_, [LeeM13]_) |br| modified version of original Kaltofen's algorithm for leading coefficient precomputation with square-free decomposition (instead of distinct variables decomposition) due to Lee is used; further adaptations are made to work in finite fields of very small cardinality; the resulting algorithm is close to [LeeM13]_, but at the same time has many differences in details:

- `MultivariateFactorization.factorInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

31. *Kaltofen's multivariate factorization Z* |____| ([Kalt85]_, [LeeM13]_) |br| (with the same modifications as for algorithm for finite fields):

- `MultivariateFactorization.factorInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

32. *Multivariate polynomial interpolation with Newton method*:

- `MultivariateInterpolation <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateInterpolation.java>`_
 

References
==========

.. [GaGe03] J von zur Gathen and J Gerhard. Modern computer algebra (2 ed.). Cambridge University Press, 2003.

.. [ShoNTL] V Shoup. NTL: A library for doing number theory. www.shoup.net/ntl

.. [GeCL92] K O Geddes, S R Czapor, G Labahn. Algorithms for Computer Algebra. 1992.

.. [Muss71] D R Musser. Algorithms for polynomial factorization, Ph.D. Thesis, University of Wisconsin, 1971.

.. [Shou95] V Shoup. A new polynomial factorization algorithm and its implementation. J. Symb. Comput., 20(4):363–397, 1995.

.. [BreK98] R P Brent and H T Kung. Fast algorithms for manipulating formal power series. J. Assoc. Comput. Math. 25:581-595, 1978

.. [Brow71] W S Brown. On Euclid’s algorithm and the computation of polynomial greatest common divisors. J. ACM, 18(4):478–504, 1971.

.. [Zipp79] R E Zippel. Probabilistic algorithms for sparse polynomials. In Proceedings of the International Symposiumon on Symbolic and Algebraic Computation, EUROSAM ’79, pages 216–226, London, UK, UK, 1979. Springer-Verlag.

.. [Zipp93] R E Zippel. Effective Polynomial Computation. Kluwer International Series in Engineering and Computer Science. Kluwer Academic Publishers, 1993.

.. [dKMW05] J de Kleine, M B Monagan, A D Wittkopf. Algorithms for the Non-monic Case of the Sparse Modular GCD Algorithm. Proceeding of ISSAC ’05, ACM Press, pp. 124-131 , 2005.

.. [Yang09] S Yang. Computing the greatest common divisor of multivariate polynomials over finite fields. Master’s thesis, Simon Fraser University, 2009.

.. [MosY73] J Moses and D Y Y Yun, "The EZGCD Algorithm," pp. 159-166 in Proc. ACM Annual Conference, (1973).

.. [Wang80] P S Wang, "The EEZ-GCD Algorithm," ACM SIGSAMBull., 14 pp. 50-60 (1980).

.. [KalM99] E Kaltofen, M. B. Monagan. On the Genericity of the Modular Polynomial GCD Algorithm. Proceeding of ISSAC ’99, ACM Press, 59-66, 1999.

.. [Bern99] L Bernardin. Factorization of Multivariate Polynomials over Finite Fields. PhD thesis, ETH Zu ̈rich, 1999.

.. [LeeM13] M M-D Lee, Factorization of multivariate polynomials,  Ph.D. thesis, University of Kaiserslautern, 2013

.. [Kalt85] E Kaltofen. Sparse Hensel lifting. In EUROCAL 85 European Conf. Comput. Algebra Proc. Vol. 2, pages 4–17, 1985.


-------

=======
License
=======

Apache License, Version 2.0 http://www.apache.org/licenses/LICENSE-2.