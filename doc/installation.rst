.. |br| raw:: html

   <br/>

.. _ref-installation:

======
Set up
======


|Rings| are currently available for Java and Scala.

To get started with Scala SBT, simply add the following dependence to your ``build.sbt`` file:

.. code-block:: scala

	libraryDependencies += "cc.redberry" % "rings.scaladsl" % "2.0"

For using |Rings| solely in Java there is Maven artifact:

.. code-block:: xml

	<dependency>
	    <groupId>cc.redberry</groupId>
	    <artifactId>rings</artifactId>
	    <version>2.0</version>
	</dependency>


|Rings| library has the following structure:

 - ``rings`` |br| the core of |Rings| library written entirely in Java. It includes:
 
 	- ``rings.bigint`` |br| arbitrary precision integers (fork of `tbuktu/bigint <https://github.com/tbuktu/bigint>`_)
 	- ``rings.primes`` |br| prime numbers including prime factorization, primality test etc.
 	- ``rings.poly.univar`` |br| univariate polynomials and algorithms with them including GCD and factorization
 	- ``rings.poly.multivar`` |br| multivariate polynomials and algorithms with them including GCD, factorization, Groebner basis etc.
 
 - ``rings.scaladsl`` |br| Scala wrappers and syntax definitions for |Rings|


Playing around
==============

To taste what |Rings| can do, one can clone the `github repo <https://github.com/PoslavskySV/rings>`_ and try |Rings| with SBT's console. Rung ``sbt console`` at the prompt:


.. code-block:: text
    
	> sbt console
	[info] Starting scala interpreter...
	[info] 
	Welcome to Scala 2.12.3 (Java HotSpot(TM) 64-Bit Server VM, Java 1.8.0_144).
	Type in expressions for evaluation. Or try :help.

    scala > import cc.redberry.rings.poly.PolynomialMethods._
    import cc.redberry.rings.poly.PolynomialMethods._

    scala> import cc.redberry.rings.scaladsl._
    import cc.redberry.rings.scaladsl._

    scala> import syntax._
    import syntax._

    scala> implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))
    ring: cc.redberry.rings.scaladsl.MultivariateRing[cc.redberry.rings.bigint.BigInteger] = Z[a, b, c]

    scala> val poly1 = ring("x + y - z").pow(8)

    scala> val poly2 = ring("x - y + z").pow(8)
    
    scala> ring show Factor(poly1 - poly2)