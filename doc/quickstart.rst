.. _ref-quickstart:

==========
Quick tour
==========


The first step in |Rings| is to import the main modules: 

.. tabs::

   .. code-tab:: scala

		import cc.redberry.rings.poly.PolynomialMethods._
		import cc.redberry.rings.scaladsl._
		import syntax._

   .. code-tab:: java

		import cc.redberry.rings.*
		import cc.redberry.rings.poly.*
		import cc.redberry.rings.poly.univar.*
		import cc.redberry.rings.poly.multivar.*

		import static cc.redberry.rings.Rings.*
		import static cc.redberry.rings.poly.PolynomialMethods.*
		


----


Univariate polynomial factorization
===================================

Univariate factorization is supported for polynomials in :math:`F[x]` where :math:`F` is either finite field or :math:`Z` or :math:`Q`.

----

Factor polynomial in :math:`Z_{17}[x]`:

.. tabs::

   .. code-tab:: scala

		// Define ring Z/17[x]
		implicit val ring = UnivariateRingZp64(17, "x")

		// parse univariate poly from string
		val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 - x^6 + 10*x^7 + x^8")

		// factorize poly
		val factors = Factor(poly)
		println(s"factorization : ${ring show factors}")


   .. code-tab:: java

		// the modulus
		long modulus = 17;
		// parse univariate poly over Z/17 from string
		UnivariatePolynomialZp64 poly = UnivariatePolynomialZp64
		    .parse("4 + 8*x + 12*x^2 + 5*x^5 - x^6 + 10*x^7 + x^8", modulus);

		// factorize poly
		FactorDecomposition<UnivariatePolynomialZp64> factors = PolynomialMethods.Factor(poly);
		System.out.println(factors);


This will give the following factorization:


.. math::
   :nowrap:

	\begin{eqnarray*}
	&& 4 + 8 x + 12 x^2 + 5 x^5 - x^6 + 10 x^7 + x^8 =  \\
	 && \qquad = (6 + x) (8 + x) (14 + x)^2 (15 + x) (7 + x + 4 x^2 + x^3) \quad \text{mod}\,\, 17
	\end{eqnarray*}


----

Coefficient rings with arbitrary large characteristic are available:

.. tabs::

   .. code-tab:: scala

		// coefficient ring Z/1237940039285380274899124357 (the next prime to 2^100)
		val cfRing = Zp(new BigInteger("1267650600228229401496703205653"))

		// The ring Z/1237940039285380274899124357[x]
		implicit val ring = UnivariateRing(cfRing, "x")

		val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 + 16*x^6 + 27*x^7 + 18*x^8")
		println(s"factorization : ${ring show Factor(poly)}")

   .. code-tab:: java

		//// coefficient ring Z/1237940039285380274899124357 (the next prime to 2^100)
		IntegersZp cfRing = Rings.Zp(new BigInteger("1267650600228229401496703205653"));

		// parse univariate poly over Z/1267650600228229401496703205653 from string
		UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial
		    .parse("4 + 8*x + 12*x^2 + 5*x^5 - x^6 + 10*x^7 + x^8", cfRing);

		// factorize poly
		FactorDecomposition<UnivariatePolynomial<BigInteger>> factors 
				= PolynomialMethods.Factor(poly);
		System.out.println(factors);


(large primes can be generated with ``BigPrimes.nextPrime(BigInteger)`` method).

This will give the following factorization:


.. math::
   :nowrap:

	\begin{eqnarray*}
	&& 4 + 8 x + 12 x^2 + 5 x^5 - x^6 + 10 x^7 + x^8 =  \\
	 && \quad = (448975734644581867134339749139 + x)\times \\
	 && \qquad \times(924109545982468663492425885021 + \\ 
	 && \qquad + 396701390689518624208584222054 x + \\ 
	 && \qquad + 671565661754860153453068172251 x^2 + x^3)\times \\ 
	 && \qquad \times(493224222589341667858050863719 + \\
	 && \qquad + 336789330550038195919829925685 x + \\ 
	 && \qquad + 636344447485090019332177467857 x^2 + \\
     && \qquad + 147109203828787380909295284273 x^3 + \\
     && \qquad \qquad + x^4) \quad \text{mod}\,\, 1267650600228229401496703205653
   \end{eqnarray*}

----


Multivariate polynomial factorization
=====================================


Factor polynomial in :math:`Z_{2}[x, y, z]`:


.. tabs::

   .. code-tab:: scala

   		// The ring Z/2[x, y, z]
		implicit val ring = MultivariateRingZp64(2, Array("x", "y", "z"))

		val (x, y, z) = ring("x", "y", "z")
		// factorize polynomial
		val factors = Factor(1 + (1 + x + y + z)**2 + (x + y + z)**4)
		println(s"factorization : ${ring show factors}")


   .. code-tab:: java

		// coefficient ring Z/2
		IntegersZp64 cfRing = new IntegersZp64(2);
		MultivariatePolynomialZp64
		        // create unit multivariate polynomial over
		        // 3 variables over Z/2 using LEX ordering
		        one = MultivariatePolynomialZp64.one(3, cfRing, MonomialOrder.LEX),
		        // create "x" polynomial
		        x = one.createMonomial(0, 1),
		        // create "y" polynomial
		        y = one.createMonomial(1, 1),
		        // create "z" polynomial
		        z = one.createMonomial(2, 1);

		// (1 + x + y + z)^2
		MultivariatePolynomialZp64 poly1 = one.copy().add(x, y, z);
		poly1 = polyPow(poly1, 2);

		// (x + y + z)^4
		MultivariatePolynomialZp64 poly2 = x.copy().add(y, z);
		poly2 = polyPow(poly2, 4);

		// 1 + (1 + x + y + z)^2 + (x + y + z)^4
		MultivariatePolynomialZp64 poly = one.copy().add(poly1, poly2);
		FactorDecomposition<MultivariatePolynomialZp64> factors = PolynomialMethods.Factor(poly);
		System.out.println(factors);


This will give the following factorization:


.. math::
   :nowrap:

	\begin{eqnarray*}
	&& 1 + (1 + x + y + z)^2 + (x + y + z)^4  = (x + y + z)^2 \, (1 + x + y + z)^2 \quad \text{mod}\,\, 2
	\end{eqnarray*}

----

Factor polynomial in :math:`Z[a, b, c]`:


.. tabs::

   .. code-tab:: scala

   		// The ring Z[a, b, c]
		implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))

		val (a, b, c) = ring("a", "b", "c")
		// factorize polynomial
		val factors = Factor(1 - (1 + a + b + c)**2 - (2 + a + b + c)**3)
		println(s"factorization : ${ring show factors}")


   .. code-tab:: java

		MultivariatePolynomial<BigInteger>
		        // create unit multivariate polynomial over
		        // 3 variables over Z using LEX ordering
		        one = MultivariatePolynomial.one(3, Rings.Z, MonomialOrder.LEX),
		        // create "a" polynomial
		        a = one.createMonomial(0, 1),
		        // create "b" polynomial
		        b = one.createMonomial(1, 1),
		        // create "c" polynomial
		        c = one.createMonomial(2, 1);

		// (1 + a + b + c)^2
		MultivariatePolynomial<BigInteger> poly1 = one.copy().add(a, b, c);
		poly1 = polyPow(poly1, 2);

		// (2 + a + b + c)**3
		MultivariatePolynomial<BigInteger> poly2 = one.copy().multiply(2).add(a, b, c);
		poly2 = polyPow(poly2, 3);

		// 1 - (1 + a + b + c)^2 - (2 + a + b + c)**3
		MultivariatePolynomial<BigInteger> poly = one.copy().subtract(poly1, poly2);
		FactorDecomposition<MultivariatePolynomial<BigInteger>> factors 
				= PolynomialMethods.Factor(poly);
		System.out.println(factors);


Will give the following factorization:


.. math::
   :nowrap:

	\begin{eqnarray*}
	&& 1 - (1 + a + b + c)^2 - (2 + a + b + c)^3 = -(1 + a + b + c) (2 + a + b + c) ( 4 + a + b + c)
	\end{eqnarray*}

----

Factor polynomial in :math:`Q[x, y, z]`:

.. tabs::

   .. code-tab:: scala

		implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))

		val poly = ring(
		  """
		    |(1/6)*y*z + (1/6)*y^3*z^2 - (1/2)*y^6*z^5 - (1/2)*y^8*z^6
		    |-(1/3)*x*z - (1/3)*x*y^2*z^2 + x*y^5*z^5 + x*y^7*z^6
		    |+(1/9)*x^2*y^2*z - (1/3)*x^2*y^7*z^5 - (2/9)*x^3*y*z
		    |+(2/3)*x^3*y^6*z^5 - (1/2)*x^6*y - (1/2)*x^6*y^3*z
		    |+x^7 + x^7*y^2*z - (1/3)*x^8*y^2 + (2/3)*x^9*y
		  """.stripMargin)

		// factorize polynomial (in this example there will be 3 factors)
		val factors = Factor(poly)
		println(s"factorization : ${ring show factors}")

   .. code-tab:: java

		MultivariatePolynomial<Rational<BigInteger>>
				poly = MultivariatePolynomial.parse(
					"(1/6)*y*z + (1/6)*y^3*z^2 - (1/2)*y^6*z^5 - (1/2)*y^8*z^6" +
			        "-(1/3)*x*z - (1/3)*x*y^2*z^2 + x*y^5*z^5 + x*y^7*z^6" +
			        "+(1/9)*x^2*y^2*z - (1/3)*x^2*y^7*z^5 - (2/9)*x^3*y*z" +
			        "+(2/3)*x^3*y^6*z^5 - (1/2)*x^6*y - (1/2)*x^6*y^3*z" +
			        "+x^7 + x^7*y^2*z - (1/3)*x^8*y^2 + (2/3)*x^9*y"
				, Rings.Q);

		System.out.println(PolynomialMethods.Factor(poly));


----

Polynomial GCD
==============

Univariate extended GCD in :math:`Z_{17}[x]`:

.. tabs::

   .. code-tab:: scala

   		// The ring Z/17[x]
		implicit var ring = UnivariateRingZp64(17, "x")

		val x = ring("x")
		
		val xgcd = PolynomialExtendedGCD(1 + x + x.pow(2) + x.pow(3), 1 + 2*x + 9*x.pow(2))
		println(s"XGCD : ${ring show xgcd}")


   .. code-tab:: java

		UnivariatePolynomialZp64
		        a = UnivariatePolynomialZ64.create(1, 1, 1, 1).modulus(17),
		        b = UnivariatePolynomialZ64.create(1, 2, 9).modulus(17);

		System.out.println(Arrays.toString(UnivariateGCD.PolynomialExtendedGCD(a, b)));


----

Multivariate GCD in :math:`Z[a, b, c]`:

.. tabs::

   .. code-tab:: scala

   		// The ring Z[a, b, c]
		implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))

		println(PolynomialGCD(
			ring("-b-b*c-b^2+a+a*c+a^2"), 
			ring("b^2+b^2*c+b^3+a*b^2+a^2+a^2*c+a^2*b+a^3")))


   .. code-tab:: java

		MultivariatePolynomial<BigInteger>
		        a = MultivariatePolynomial.parse("-b-b*c-b^2+a+a*c+a^2", Rings.Z),
		        b = MultivariatePolynomial.parse("b^2+b^2*c+b^3+a*b^2+a^2+a^2*c+a^2*b+a^3", Rings.Z);

		System.out.println(PolynomialMethods.PolynomialGCD(a, b));


----

Constructing arbitrary rings
============================

Polynomial rings over :math:`Z` and :math:`Q`:

.. tabs::

	.. code-tab:: scala

		// Ring Z[x]
		UnivariateRing(Z, "x")
		// Ring Z[x, y, z]
		MultivariateRing(Z, Array("x", "y", "z"))
		// Ring Q[a, b, c]
		MultivariateRing(Q, Array("a", "b", "c"))

 	.. code-tab:: java

		// Ring Z[a]
		Rings.UnivariateRing(Rings.Z);
		// Ring Z[a, b, c]
		Rings.MultivariateRing(3, Rings.Z);
		// Ring Q[a, b, c]
		Rings.MultivariateRing(3, Rings.Q);


Polynomial rings over :math:`Z_p`:

.. tabs::

	.. code-tab:: scala

		// Ring Z/3[x] (64 indicates that machine numbers are used in the basis)
		UnivariateRingZp64(3, "x")
		// Ring Z/3[x, y, z]
		MultivariateRingZp64(3, Array("x", "y", "z"))
		// Ring Z/p[x, y, z] with p = 2^107 - 1 (Mersenne prime)
		MultivariateRing(Zp(BigInt(2).pow(107) - 1), Array("x", "y", "z"))

	.. code-tab:: java

		// Ring Z/3[a] (64 indicates that machine numbers are used in the basis)
		Rings.UnivariateRingZp64(3);
		// Ring Z/3[a, b, c]
		Rings.MultivariateRingZp64(3, 3);
		// Ring Z/p[a, b, c] with p = 2^107 - 1 (Mersenne prime)
		Rings.MultivariateRing(3, Rings.Zp(BigInteger.ONE.shiftLeft(107).decrement()));


Galois fields:

.. tabs::

   .. code-tab:: scala

		// Galois field with cardinality 7^10 
		// (irreducible polynomial will be generated automatically)
		GF(7, 10, "x")
		// GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
		GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")

   .. code-tab:: java

		// Galois field with cardinality 7^10 
		// (irreducible polynomial will be generated automatically)
		Rings.GF(7, 10);
		// GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
		Rings.GF(UnivariatePolynomialZ64.create(1, 3, 1, 2).modulus(7));


Fractional fields:

.. tabs::

   .. code-tab:: scala

		// Field of fractions of univariate polynomials Z[x]
		Rationals(UnivariateRing(Z, "x"))
		// Field of fractions of multivariate polynomials Z/19[x, y, z]
		Rationals(MultivariateRingZp64(19, Array("x", "y", "z")))

   .. code-tab:: java

		// Field of fractions of univariate polynomials Z[a]
		Rings.Rationals(Rings.UnivariateRing(Rings.Z));
		// Field of fractions of multivariate polynomials Z/19[a, b, c]
		Rings.Rationals(Rings.MultivariateRingZp64(3, 19));


----

Ring of univariate polynomials over elements of Galois field :math:`GF(7^{3})[x]`:

.. tabs::

   .. code-tab:: scala

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


   .. code-tab:: java

		// Elements of GF(7^3) are represented as polynomials
		// modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
		FiniteField<UnivariatePolynomialZp64> cfRing 
		        = Rings.GF(UnivariatePolynomialZ64.create(1, 3, 1, 2).modulus(7));
		assert cfRing.characteristic().intValue() == 7;
		assert cfRing.cardinality().intValue() == 343;

		// Ring GF(7^3)[a]
		UnivariateRing<UnivariatePolynomial<UnivariatePolynomialZp64>>
		        ring = Rings.UnivariateRing(cfRing);

		// Coefficients of polynomials in GF(7^3)[a] are elements of GF(7^3)
		UnivariatePolynomial<UnivariatePolynomialZp64> 
		        poly = ring.parse("1 - (1 - z^3) * x^6 + (1 - 2*z) * x^33 + x^66");

		// factorize poly (in this examples there will be 9 factors)
		FactorDecomposition<UnivariatePolynomial<UnivariatePolynomialZp64>> factors 
		        = PolynomialMethods.Factor(poly);
		System.out.println(factors);


Ring of multivariate polynomials over elements of Galois field :math:`GF_{7^{3}}[x, y, z]`:

.. tabs::

   .. code-tab:: scala

		// Elements of GF(7^3) are represented as polynomials
		// over "z" modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
		val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
		// Ring GF(7^3)[x]
		implicit val ring = MultivariateRing(cfRing, Array("a", "b", "c"))

		// Coefficients of polynomials in GF(7^3)[x] are elements of GF343
		val poly = ring("1 - (1 - z^3) * a^6*b + (1 - 2*z) * c^33 + a^66")


   .. code-tab:: java

		// Elements of GF(7^3) are represented as polynomials
		// modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
		FiniteField<UnivariatePolynomialZp64> cfRing
		        = Rings.GF(UnivariatePolynomialZ64.create(1, 3, 1, 2).modulus(7));
		assert cfRing.characteristic().intValue() == 7;
		assert cfRing.cardinality().intValue() == 343;

		// Ring GF(7^3)[a, b, c]
		MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>>
		        ring = Rings.MultivariateRing(3, cfRing);

		// Coefficients of polynomials in GF(7^3)[a, b, c] are elements of GF(7^3)
		MultivariatePolynomial<UnivariatePolynomialZp64>
		        poly = ring.parse("1 - (1 - z^3) * a^6*b + (1 - 2*z) * c^33 + a^66");