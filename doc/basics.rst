.. _ref-basicconcepts:

==============
Basic concepts
==============



Integers
========

There are two basic types of integer numbers that we have to deal with when doing algebra in computer: machine integers and arbitrary-precision integers. For the machine integers the Java's primitive 64-bit ``long`` type is used (since most modern CPUs are 64-bit). Internally |Rings| use machine numbers for representing integers modulo prime numbers less than :math:`2^{64}` which is done for performance reasons (see :ref:`ref-machine-arithmetic`). For the arbitrary-precision integers |Rings| use improved ``BigInteger`` class `github.com/tbuktu/bigint <https://github.com/tbuktu/bigint>`_ (`cc.redberry.rings.bigint.BigInteger`_) instead of built-in ``java.math.BigInteger``. The improved ``BigInteger`` has Schönhage-Strassen multiplication and Barrett division algorithms for large integers which is a significant performance improvement in comparison to native Java's implementation.


.. _cc.redberry.rings.bigint.BigInteger: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/bigint/BigInteger.java

Prime numbers
"""""""""""""

In many applications it is necessary to test primality of integer number (``isPrime(number)``) or to generate some prime numbers at random (``nextPrime(number)``). This is realized in the following two classes:

 - `SmallPrimes`_ for numbers less than :math:`2^{32}`. It uses Miller-Rabin probabilistic primality test for int type in such a way that a result is always guaranteed (code is adapted from `Apache Commons Math <http://commons.apache.org/proper/commons-math/>`_).
 - `BigPrimes`_ for arbitrary large numbers. It switches between Pollard-Rho, Pollard-P1 and Quadratic Sieve algorithms for prime factorization and also uses probabilistic Miller-Rabin test and strong Lucas test for primality testing.


The following examples give some illustrations:


.. tabs::

   .. code-tab:: java

		int intNumber = 1234567;
		// false
		boolean primeQ = SmallPrimes.isPrime(intNumber);
		// 1234577
		int intPrime = SmallPrimes.nextPrime(intNumber);
		// [127, 9721]
		int[] intFactors = SmallPrimes.primeFactors(intNumber);


		long longNumber = 12345671234567123L;
		// false
		primeQ = BigPrimes.isPrime(longNumber);
		// 12345671234567149
		long longPrime = BigPrimes.nextPrime(longNumber);
		// [1323599, 9327350077]
		long[] longFactors = BigPrimes.primeFactors(longNumber);


		BigInteger bigNumber = new BigInteger("321536584276145124691487234561832756183746531874567");
		// false
		primeQ = BigPrimes.isPrime(bigNumber);
		// 321536584276145124691487234561832756183746531874827
		BigInteger bigPrime = BigPrimes.nextPrime(bigNumber);
		// [3, 29, 191, 797359, 1579057, 14916359, 1030298906727233717673336103]
		List<BigInteger> bigFactors = BigPrimes.primeFactors(bigNumber);


.. _SmallPrimes: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/primes/SmallPrimes.java
.. _BigPrimes: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/primes/BigPrimes.java


.. _ref-machine-arithmetic:

Modular arithmetic with machine integers
========================================

There is one special ring --- ring :math:`Z_p` of integers modulo prime number :math:`p < 2^{64}` --- which is used in the basis of many fundamental algorithms. In contrast to :math:`Z_p` with arbitrary large characteristic, for characteristic that fits into 64-bit word one can use machine integers to significantly speed up basic math operations. Operations in :math:`Z_p` require applying ``mod`` operation which in turn implies integer division. Integer division is a very slow CPU instruction; and what is more important is that it breaks CPU pipelining. Hopefully, operations in :math:`Z_p` imply taking ``mod`` with a fixed modulus :math:`p`, which allows to make some precomputation beforehand and then reduce integer divisions to multiplications that are over a magnitude times faster. The details of this trick can be found in `Hacker's Delight <http://www.hackersdelight.org>`_. |Rings| use `libdivide4j`_ library for fast integer division with precomputation which is ported from the well known C/C++ `libdivide`_ library. With this precomputation the ``mod`` operation becomes several times faster than the native CPU instruction, which boosts the overall performance of many of |Rings| algorithms in more than 3 times.

.. _libdivide4j: https://github.com/PoslavskySV/libdivide4j/

.. _libdivide: https://libdivide.com

The ring :math:`Z_p` with :math:`p < 2^{64}` is implemented in `IntegersZp64`_ class (while `IntegersZp`_ implements :math:`Z_p` with arbitrary large characteristic). `IntegersZp64`_ defines all arithmetic operations in :math:`Z_p`:

.. tabs::

   .. code-tab:: java

		// Z/p with p = 2^7 - 1 (Mersenne prime)
		IntegersZp64 field = new IntegersZp64(127);
		//     1000 = 111 mod 127
		assert field.modulus(1000) == 111;
		// 100 + 100 = 73 mod 127
		assert field.add(100, 100) == 73;
		//  12 - 100 = 39 mod 127
		assert field.subtract(12, 100) == 39;
		//  55 * 78  = 73 mod 127
		assert field.multiply(55, 78) == 99;
		//   1 / 43  = 65 mod 127
		assert field.reciprocal(43) == 65;

It is worst to mention, that multiplication defined in `IntegersZp64`_ is especially fast when characteristic is less than :math:`2^{32}`: in this case multiplication of two numbers fits the machine 64-bit word, while in the opposite case Montgomery reduction will be used.


.. tabs::

   .. code-tab:: java

   		// Z/p with p = 2^31 - 1 (Mersenne prime) - fits 32-bit word
		IntegersZp64 field32 = new IntegersZp64((1L << 31) - 1L);
		// does cause long overflow - fast 
		assert field32.multiply(0xabcdef12, 0x12345678) == 0x7e86a4d6;


		// Z/p with p = 2^61 - 1 (Mersenne prime) - doesn't fit 32-bit word
		IntegersZp64 field64 = new IntegersZp64((1L << 61) - 1L);
		// cause long overflow - Montgomery reduction will be used - no so fast 
		assert field64.multiply(0x0bcdef1234567890L, 0x0234567890abcdefL) == 0xf667077306fd7a8L;



**Implementation note:** unfortunately, the price that we pay for fast arithmetic with machine integers is that `IntegersZp64`_ stands separately from the elegant type hierarchy of generic rings implemented in |Rings| (see section :ref:`ref-rings`); that is because Java doesn't support generics with primitive types. This leads to that some of the fundamental algorithms have two implementations -- one for rings over generic elements and one for `IntegersZp64`_ over ``longs``.


.. _IntegersZp64: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/IntegersZp64.java
.. _IntegersZp: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/IntegersZp.java


.. _ref-rings:

Rings
=====

The concept of mathematical ring is implemented in the generic interface `Ring<E>`_ which defines all basic algebraic operations over the elements of type ``E``. The simplest example is the ring of integers :math:`Z` (`Rings.Z`_), which operates with ``BigInteger`` instances and simply delegates all operations like ``+`` or ``*`` to methods of class ``BigInteger``. A little bit more complicated ring is a ring of integers modulo some number :math:`Z_p`:

.. tabs::

   .. code-tab:: java

		// The ring Z/17
		Ring<BigInteger> ring = Rings.Zp(BigInteger.valueOf(17));
		
		//     103 = 1 mod 17 
		BigInteger el  = ring.valueOf(BigInteger.valueOf(103));
		assert  el.intValue() == 1;
		
		// 99 + 88 = 0 mod 17
		BigInteger add = ring.add(BigInteger.valueOf(99),
		                          BigInteger.valueOf(88));
		assert add.intValue() == 0;

		// 99 * 77 = 7 mod 17
		BigInteger mul = ring.multiply(BigInteger.valueOf(99),
		                               BigInteger.valueOf(77));
		assert mul.intValue() == 7;

		// 1  / 99 = 11 mod 17
		BigInteger inv = ring.reciprocal(BigInteger.valueOf(99));
		assert inv.intValue() == 11;


In fact the interface `Ring<E>`_ defines algebraic operations inherent both for *GCD domains*, *Euclidean rings* and *Fields*. These operations can be summarized in the following methods from `Ring<E>`_:


.. tabs::

   .. code-tab:: java

		// Methods from Ring<E> interface:

		// GCD domain operation:
		E gcd(E a, E b);

		// Euclidean ring operation:
		E[] divideAndRemainder(E dividend, E divider);

		// Field operation:
		E reciprocal(E element);

In the case when a particular ring is (e.g. :math:`Z`) is not a field, the invocation of corresponding method (``reciprocal``) will produce ``ArithmeticException``. Each `Ring<E>`_ implementation provides the information about its mathematical origin and all properties like cardinality, characteristic etc. Additionally it defines ``parse(String)`` method to convert strings into ring elements:


.. tabs::

   .. code-tab:: java

		// Z is not a field
		assert Rings.Z.isEuclideanRing();
		assert !Rings.Z.isField();
		assert !Rings.Z.isFinite();

		// Q is an infinite field
		assert Rings.Q.isField();
		assert !Rings.Q.isFinite();
		assert Rings.Q.parse("2/3").equals(
			new Rational<>(Rings.Z, BigInteger.valueOf(2), BigInteger.valueOf(3)));

		// GF(2^10) is a finite field
		FiniteField<UnivariatePolynomialZp64> gf = Rings.GF(2, 10);
		assert gf.isField();
		assert gf.isFinite();
		assert gf.characteristic().intValue() == 2;
		assert gf.cardinality().intValue() == 1 << 10;
		System.out.println(gf.parse("1 + z + z^10"));

		// Z/3[x] is Euclidean ring but not a field
		UnivariateRing<UnivariatePolynomialZp64> zp3x = Rings.UnivariateRingZp64(3);
		assert zp3x.isEuclideanRing();
		assert !zp3x.isField();
		assert !zp3x.isFinite();
		assert zp3x.characteristic().intValue() == 3;
		assert zp3x.parse("1 + 14*x + 15*x^10").equals(
			UnivariatePolynomialZ64.create(1, 2).modulus(3));


Examples of rings
"""""""""""""""""

The shortcut methods for different rings are placed in `cc.redberry.rings.Rings`_ class. Below is the list of basic rings defined in |Rings|:

+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| Ring                                   | Description                                                         | Code in Rings                         |
+========================================+=====================================================================+=======================================+
| :math:`Z`                              | Ring of integers                                                    | ``Z``                                 |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Q`                              | Field of rationals                                                  | ``Q``                                 |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Z_p`                            | Integers modulo :math:`p`                                           | ``Zp(p)``                             |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Z_p` with :math:`p < 2^{64}`    | Integers modulo :math:`p < 2^{64}`                                  | ``Zp64(p)`` [*]_                      |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`GF(p^q)`                        | Galois field with cardinality :math:`p^q`                           | ``GF(p, q)`` or ``GF(irred)``         |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Frac(R)`                        | Field of fractions of an integral domain :math:`R`                  | ``Frac(R)``                           |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`R[x]`                           | Univariate polynomial ring over                                     | ``UnivariatePolynomials(R)``          |
|                                        | coefficient ring :math:`R`                                          |                                       |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Z_p[x]` with :math:`p < 2^{64}` | Univariate polynomial ring over                                     | ``UnivariatePolynomialsZp64(p)``      |
|                                        | coefficient ring :math:`Z_p` with :math:`p < 2^{64}`                |                                       |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`R[x_1, \dots, x_N]`             | Multivariate polynomial ring with exactly :math:`N`                 | ``MultivariatePolynomials(N, R)``     |
|                                        | variables over coefficient ring :math:`R`                           |                                       |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+
| :math:`Z_p[x_1, \dots, x_N]`           | Multivariate polynomial ring with exactly :math:`N`                 | ``MultivariatePolynomialsZp64(N, p)`` |
| with :math:`p < 2^{64}`                | variables over coefficient ring :math:`Z_p` with :math:`p < 2^{64}` |                                       |
+----------------------------------------+---------------------------------------------------------------------+---------------------------------------+


.. [*] Class `IntegersZp64`_ which represents :math:`Z_p` with :math:`p < 2^{64}` does not inherit `Ring<E>`_ interface (see :ref:`ref-machine-arithmetic`)


Scala DSL defines a wrapper class `Ring[E]`_ with implicit conversion to Java's `Ring<E>`_ for the reasons described below (see :ref:`ref-basics-polynomials`). So in Scala constructor methods from `cc.redberry.rings.Rings`_ are defined in `cc.redberry.rings.scaladsl.Rings`_.


Galois fields
^^^^^^^^^^^^^

Galois field :math:`GF(p^q)` with prime characteristic :math:`p` and cardinality :math:`p^q` can be can be created by specifying :math:`p` and :math:`q` in which case the irreducible polynomial will be generated automatically or by explicitly specifying the irreducible:

.. tabs::

   .. code-tab:: scala

		// Galois field GF(7^10) represented by univariate polynomials 
		// in variable "z" over Z/7 modulo some irreducible polynomial
		// (irreducible polynomial will be generated automatically)
		GF(7, 10, "z")
		// GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
		GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")

   .. code-tab:: java

		// Galois field GF(7^10)
		// (irreducible polynomial will be generated automatically)
		GF(7, 10);
		// GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
		GF(UnivariatePolynomialZ64.create(1, 3, 1, 2).modulus(7));

Galois fields with arbitrary large characteristic are available:

.. tabs::

	.. code-tab:: scala

		// Mersenne prime 2^107 - 1
		val characteristic : BigInteger = BigInt(2).pow(107) - 1
		// Galois field GF((2^107 - 1) ^ 16)
		implicit val field = GF(characteristic, 16, "z")
		
		assert(field.cardinality() == characteristic.pow(16))
		

	.. code-tab:: java

		// Mersenne prime 2^107 - 1
		BigInteger characteristic = BigInteger.ONE.shiftLeft(107).decrement();
		// Galois field GF((2^107 - 1) ^ 16)
		FiniteField<UnivariatePolynomial<BigInteger>> field = GF(characteristic, 16);

		assert(field.cardinality().equals(characteristic.pow(16)));


Implementation of Galois fields uses precomputed inverses for fast division with Newton iterations (see ``fastDivisionPreConditioning`` in `UnivariateDivision`_) which allows to achieve assymptotically fast performance.


Fields of fractions
^^^^^^^^^^^^^^^^^^^

Field of fractions can be defined over any integral domain :math:`R`. The simplest example is the field :math:`Q` of fractions over :math:`Z`:

.. tabs::

	.. code-tab:: scala

		implicit val field = Frac(Z) // the same as Q

		assert( field("13/6") == field("2/3") + field("3/2") ) 
		

	.. code-tab:: java

		Rationals<BigInteger> field = Frac(Z); // the same as Q

		Rational<BigInteger> a = field.parse("13/6");
		Rational<BigInteger> b = field.parse("2/3");
		Rational<BigInteger> c = field.parse("3/2");

		assert a.equals(field.add(b, c));


The common GCD is automatically canceled in the numerator and denominator. Fractions may be defined over any GCD ring. For example, :math:`Frac(Z[x, y, z])` -- rational functions over :math:`x`, :math:`y` and :math:`z`:


.. tabs::

	.. code-tab:: scala

		val ring = MultivariateRing(Z, Array("x", "y", "z"))
		implicit val field = Frac(ring)

		val a = field("(x + y + z)/(1 - x - y)")
		val b = field("(x^2 - y^2 + z^2)/(1 - x^2 - 2*x*y - y^2)")

		println(a + b)		

	.. code-tab:: java

		Ring<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(3, Z);
		Ring<Rational<MultivariatePolynomial<BigInteger>>> field = Frac(ring);

		Rational<MultivariatePolynomial<BigInteger>> 
				a = field.parse("(x + y + z)/(1 - x - y)"),
				b = field.parse("(x^2 - y^2 + z^2)/(1 - x^2 - 2*x*y - y^2)");

		System.out.println(field.add(a, b));


.. _Ring<E>: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/Ring.java

.. _Ring[E]: https://github.com/PoslavskySV/rings/blob/develop/rings.scaladsl/src/main/scala/cc/redberry/rings/scaladsl/Rings.scala

.. _Rings.Z: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/Rings.java#L30

.. _cc.redberry.rings.Rings: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/Rings.java

.. _cc.redberry.rings.scaladsl.Rings: https://github.com/PoslavskySV/rings/blob/develop/rings.scaladsl/src/main/scala/cc/redberry/rings/scaladsl/Rings.scala

.. _UnivariateDivision: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateDivision.java


.. _ref-basics-polynomials:

Polynomials and polynomial rings
================================


|Rings| have separate implementation of univariate (dense) and multivariate (sparse) polynomials. Polynomials over :math:`Z_p` with :math:`p < 2^{64}` are also implemented separately and specifically optimized (coefficients are represented as primitive machine integers instead of generic templatized objects and fast modular arithmetic is used, see :ref:`ref-machine-arithmetic`). Below the type hierarchy of polynomial classes is shown:

.. figure:: _static/PolyUML.png
   :scale: 100%
   :align: center




.. _IPolynomial<PolyType>: https://github.com/PoslavskySV/rings/blob/develop/rings/src/main/java/cc/redberry/rings/poly/IPolynomial.java

Scala DSL
=========