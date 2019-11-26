[![image](https://img.shields.io/circleci/project/github/PoslavskySV/rings.svg?style=flat)](https://circleci.com/gh/PoslavskySV/rings)
[![image](https://readthedocs.org/projects/rings/badge/?version=latest)](https://rings.readthedocs.io)
[![image](http://www.javadoc.io/badge/cc.redberry/rings.svg)](http://www.javadoc.io/doc/cc.redberry/rings)
[![image](http://www.javadoc.io/badge/cc.redberry/rings.scaladsl_2.12.svg?label=scaladoc)](http://www.javadoc.io/doc/cc.redberry/rings.scaladsl_2.12)
[![image](https://img.shields.io/maven-central/v/cc.redberry/rings/2.svg?style=flat)](https://search.maven.org/#artifactdetails%7Ccc.redberry%7Crings%7C2.5.5%7Cjar)
[![image](https://img.shields.io/maven-central/v/cc.redberry/rings.scaladsl_2.12/2.svg?style=flat)](https://search.maven.org/#artifactdetails%7Ccc.redberry%7Crings.scaladsl_2.12%7C2.5.5%7Cjar)
[![image](https://img.shields.io/badge/License-Apache%202.0-blue.svg?style=flat)](https://opensource.org/licenses/Apache-2.0)

Rings: efficient Java/Scala library for polynomial rings
========================================================

Rings is an efficient lightweight library for commutative algebra. Polynomial arithmetic, GCDs, polynomial factorization and Groebner bases are implemented with the use of modern asymptotically fast algorithms. Rings can be easily interacted or embedded in applications via simple API with fully typed hierarchy of algebraic structures and algorithms for commutative algebra. As well, an interactive REPL is also provided. The use of Scala language brings a quite novel powerful strongly typed functional programming model allowing to write short, expressive and fast code for applications. At the same time Rings shows one of the best or even unmatched in some cases performance among existing software for algebraic calculations.

The key features of Rings include:

> -  [Rings →](http://rings.readthedocs.io/en/latest/guide.html#ref-rings) Integers, fractions, finite and algebraic fields, multiple field extensions, polynomial rings and more
> -  [Polynomials →](http://rings.readthedocs.io/en/latest/guide.html#ref-basics-polynomials) Efficient univariate and multivariate polynomials over arbitrary coefficient rings
> -  [Polynomial GCD →](http://rings.readthedocs.io/en/latest/guide.html#ref-polynomial-methods) Highly performant polynomial GCD over arbitrary coefficient domains
> -  [Univariate and multivariate polynomial factorization →](http://rings.readthedocs.io/en/latest/guide.html#ref-multivariate-factorization) Highly performant polynomial factorization over almost arbitrary rings
> -  [Ideals and Gröbner bases →](http://rings.readthedocs.io/en/latest/guide.html#ref-ideals) Polynomial ideals and efficient algorithms for Gröbner bases
> -  [Scala DSL →](http://rings.readthedocs.io/en/latest/guide.html#ref-scala-dsl) Powerful domain specific language in Scala
> -  [Fast →](https://github.com/PoslavskySV/rings.benchmarks) Really fast library suitable for real-world computational challenges

The full documentation is available at [<http://rings.readthedocs.io>](https://rings.readthedocs.io).

A more academic description of the library can be found in:

> Stanislav Poslavsky, _Rings: An efficient Java/Scala library for polynomial rings_, Computer Physics Communications, Volume 235, 2019, Pages 400-413, [doi:10.1016/j.cpc.2018.09.005](https://doi.org/10.1016/j.cpc.2018.09.005)

(please, cite this paper if you use Rings)


Set up
------

### Interactive Rings shell and Rings scripts

To taste what Rings can do, one can try interactive session with [Ammonite REPL](http://ammonite.io). You can install Rings<i>.repl</i> with Homebrew:

``` bash
$ brew install PoslavskySV/rings/rings.repl
```

or just by typing the following commands at the prompt:

``` bash
$ sudo sh -c '(echo "#!/usr/bin/env sh" && curl -L https://github.com/lihaoyi/Ammonite/releases/download/1.1.2/2.12-1.1.2) > /usr/local/bin/amm && chmod +x /usr/local/bin/amm'
$ sudo sh -c 'curl -L -o /usr/local/bin/rings.repl https://git.io/vd7EY && chmod +x /usr/local/bin/rings.repl'
```

Now run Rings<i>.repl</i>:

``` scala
$ rings.repl
Loading...
Rings 2.5.5: efficient Java/Scala library for polynomial rings

@ implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))
ring: MultivariateRing[IntZ] = MultivariateRing(Z, Array("x", "y", "z"), LEX)

@ val poly1 = ring("x + y - z").pow(8) 
poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*...

@ val poly2 = ring("x - y + z").pow(8) 
poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*...

@ Factor(poly1 - poly2)
res13: FactorDecomposition[MultivariatePolynomial[IntZ]] = 
       16*(x)*((-1)*z+y)
       *(z^4-4*y*z^3+6*y^2*z^2-4*y^3*z+y^4+6*x^2*z^2-12*x^2*y*z+6*x^2*y^2+x^4)
       *(z^2-2*y*z+y^2+x^2)
```

Additionally, Rings<i>.repl</i> can be used to run scripts with Rings code:

```
$ rings.repl myRingsScript.sc
```

### Java/Scala library

Rings is currently available for Java and Scala. To get started with Scala SBT, simply add the following dependence to your `build.sbt` file:

``` scala
libraryDependencies += "cc.redberry" %% "rings.scaladsl" % "2.5.5"
```

For using Rings solely in Java there is Maven artifact:

``` scala
<dependency>
    <groupId>cc.redberry</groupId>
    <artifactId>rings</artifactId>
    <version>2.5.5</version>
</dependency>
```

### Development version

Download latest Rings from the develop:
```bash
git clone https://github.com/PoslavskySV/rings.git
cd rings
```

Install Java artifact locally:
```bash
cd rings
mvn install -DskipTests
```

Install Rings.scaladsl locally:
```bash
cd rings.scaladsl
sbt publishLocal
```

To run a simple REPL run e.g.:
```bash
cd rings.scaladsl
sbt console
```

```scala
@
import cc.redberry.rings
import cc.redberry.rings.primes.{SmallPrimes, BigPrimes}
import rings.{bigint, primes, linear, poly}
import poly.{univar, multivar}
import poly.PolynomialMethods._
import multivar.MonomialOrder._
import multivar.GroebnerMethods
import rings.scaladsl._
import util._
import syntax._

@ implicit val ring = MultivariateRing(Z, Array("x", "y", "z"))
ring: MultivariateRing[IntZ] = MultivariateRing(Z, Array("x", "y", "z"), LEX)

@ val poly1 = ring("x + y - z").pow(8) 
poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*...

@ val poly2 = ring("x - y + z").pow(8) 
poly1: MultivariatePolynomial[IntZ] = z^8-8*y*z^7+28*y^2*z^6-56*y^3*z^5+70*...

@ Factor(poly1 - poly2)
res13: FactorDecomposition[MultivariatePolynomial[IntZ]] = 
       16*(x)*((-1)*z+y)
       *(z^4-4*y*z^3+6*y^2*z^2-4*y^3*z+y^4+6*x^2*z^2-12*x^2*y*z+6*x^2*y^2+x^4)
       *(z^2-2*y*z+y^2+x^2)
```



Examples: rings, ideals, Gröbner bases, GCDs & factorization
-------------------------------------------------------------

Below examples can be evaluated directly in the Rings<i>.repl</i>. If using Rings in Scala, the following preambula will import all required things from Rings library:

``` scala
import cc.redberry.rings

import rings.poly.PolynomialMethods._
import rings.scaladsl._
import syntax._
```

Java examples can be found in the [complete documentation pages](https://rings.readthedocs.io).

------------------------------------------------------------------------

### Some built-in rings


Polynomial rings over *Z* and *Q*:

``` scala
// Ring Z[x]
UnivariateRing(Z, "x")
// Ring Z[x, y, z]
MultivariateRing(Z, Array("x", "y", "z"))
// Ring Q[a, b, c]
MultivariateRing(Q, Array("a", "b", "c"))
```

Polynomial rings over *Z_p*:

``` scala
// Ring Z/3[x]
UnivariateRingZp64(3, "x")
// Ring Z/3[x, y, z]
MultivariateRingZp64(3, Array("x", "y", "z"))
// Ring Z/p[x, y, z] with p = 2^107 - 1 (Mersenne prime)
MultivariateRing(Zp(Z(2).pow(107) - 1), Array("x", "y", "z"))
```

Galois fields:

``` scala
// Galois field with cardinality 7^10 
// (irreducible polynomial will be generated automatically)
GF(7, 10, "x")
// GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
```

Fields of rational functions:

``` scala
// Field of fractions of univariate polynomials Z[x]
Frac(UnivariateRing(Z, "x"))
// Field of fractions of multivariate polynomials Z/19[x, y, z]
Frac(MultivariateRingZp64(19, Array("x", "y", "z")))
```

### Univariate polynomials

Some algebra in Galois field *GF(17,9)*:

``` scala
// Galois field GF(17, 9) with irreducible 
// poly in Z/17[t] generated automaticaly
implicit val ring = GF(17, 9, "t")

// pick some random field element
val a = ring.randomElement()
// raise field element to the power of 1000
val b = a.pow(1000)
// reciprocal of field element
val c = 1 / b

assert ( b * c === 1)

// explicitly parse field element from string:
// input poly will be automatically converted to
// element of GF(17, 9) (reduced modulo field generator)
val d = ring("1 + t + t^2 + t^3 + 15 * t^999")
// do some arbitrary math ops in the field
val some = a / (b + c) + a.pow(6) - a * b * c * d
```

------------------------------------------------------------------------

Extended GCD in *Z_{17}[x]*:

``` scala
// polynomial ring Z/17[x]
implicit val ring = UnivariateRingZp64(17, "x")
// parse ring element
val x = ring("x")

// construct some polynomials
val poly1 = 1 + x + x.pow(2) + x.pow(3)
val poly2 = 1 + 2 * x + 9 * x.pow(2)

// compute (gcd, s, t) such that s * poly1 + t * poly2 = gcd
val Array(gcd, s, t) = PolynomialExtendedGCD(poly1, poly2)
assert (s * poly1 + t * poly2 == gcd)

println((gcd, s, t))
```

------------------------------------------------------------------------

Factor polynomial in *Z_{17}[x]*:

``` scala
// polynomial ring Z/17[x]
implicit val ring = UnivariateRingZp64(17, "x")x

// parse polynomial from string
val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 - x^6 + 10*x^7 + x^8")

// factorize poly
val factors = Factor(poly)

println(factors)
```

Coefficient rings with arbitrary large characteristic are available:


``` scala
// coefficient ring Z/1237940039285380274899124357 (the next prime to 2^100)
val modulus = Z("1267650600228229401496703205653")
val cfRing  = Zp(modulus)

// ring Z/1237940039285380274899124357[x]
implicit val ring = UnivariateRing(cfRing, "x")
val poly = ring("4 + 8*x + 12*x^2 + 5*x^5 + 16*x^6 + 27*x^7 + 18*x^8")

// factorize poly
println(Factor(poly))
```

(large primes can be generated with `BigPrimes.nextPrime` method, see [Prime numbers](http://rings.readthedocs.io/en/latest/guide.html#ref-primes)).


------------------------------------------------------------------------

Ring of univariate polynomials over elements of Galois field *GF(7,3)[x]*:

``` scala
// elements of coefficient field GF(7,3) are represented as polynomials
// over "z" modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")

assert(cfRing.characteristic().intValue() == 7)
assert(cfRing.cardinality().intValue() == 343)

// polynomial ring GF(7^3)[x]
implicit val ring = UnivariateRing(cfRing, "x")

// parse poly in GF(7^3)[x] from string
// coefficients of polynomials in GF(7,3)[x] are elements
// of GF(7,3) that is polynomials over "z"
val poly = ring("1 - (1 - z^3) * x^6 + (1 - 2*z) * x^33 + x^66")

// factorize poly
val factors = Factor(poly)
println(s"${ring show factors}")
```

### Multivariate polynomials

Some math with multivariate polynomials from *Z[x, y, z]*:

``` scala
// ring Z[x, y, z]
implicit val ring = MultivariateRing(Z, Array("x", "y", "z")) 
// parse some ring elements
val (x, y, z) = ring("x", "y", "z") 

// construct some polynomials using different math ops
val a = (x + y + z).pow(2) - 1 
val b = (x - y - z - 1).pow(2) + x + y + z - 1 
val c = (a + b + 1).pow(9) - a - b - 1

// reduce c modulo a and b (multivariate division with remainder)
val (div1, div2, rem) = c /%/% (a, b)
```

------------------------------------------------------------------------

Multivariate GCD in *Z[a, b, c]*:

``` scala
// ring Z[a, b, c]
implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))

// parse polynomials from strings
val poly1 = ring("-b-b*c-b^2+a+a*c+a^2")
val poly2 = ring("b^2+b^2*c+b^3+a*b^2+a^2+a^2*c+a^2*b+a^3")

// compute multivariate GCD
val gcd   = PolynomialGCD(poly1, poly2)
assert (poly1 % gcd === 0)
assert (poly2 % gcd === 0)
println(gcd)
```

------------------------------------------------------------------------

Factor polynomial in *Z_{2}[x, y, z]*:

``` scala
// ring Z/2[x, y, z]
implicit val ring = MultivariateRingZp64(2, Array("x", "y", "z"))
val (x, y, z) = ring("x", "y", "z")

// factorize poly
val factors = Factor(1 + (1 + x + y + z).pow(2) + (x + y + z).pow(4))
println(factors)
```

------------------------------------------------------------------------

Factor polynomial in *Z[a, b, c]*:

``` scala
// ring Z[a, b, c]
implicit val ring = MultivariateRing(Z, Array("a", "b", "c"))
val (a, b, c) = ring("a", "b", "c")

// factorize poly
val factors = Factor(1 - (1 + a + b + c).pow(2) - (2 + a + b + c).pow(3))
println(ring show factors)
```

------------------------------------------------------------------------

Factor polynomial in *Q[x, y, z]*:

``` scala
// ring Q[x, y, z]
implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))

// parse some poly from string
val poly = ring(
  """
    |(1/6)*y*z + (1/6)*y^3*z^2 - (1/2)*y^6*z^5 - (1/2)*y^8*z^6
    |-(1/3)*x*z - (1/3)*x*y^2*z^2 + x*y^5*z^5 + x*y^7*z^6
    |+(1/9)*x^2*y^2*z - (1/3)*x^2*y^7*z^5 - (2/9)*x^3*y*z
    |+(2/3)*x^3*y^6*z^5 - (1/2)*x^6*y - (1/2)*x^6*y^3*z
    |+x^7 + x^7*y^2*z - (1/3)*x^8*y^2 + (2/3)*x^9*y
  """.stripMargin)

// factorize poly
val factors = Factor(poly)
println(factors)
```

------------------------------------------------------------------------

Ring of multivariate polynomials over elements of Galois field *GF(7,3)[x, y, z]*:

``` scala
// elements of GF(7,3) are represented as polynomials
// over "z" modulo irreducible polynomial "1 + 3*z + z^2 + z^3"
val cfRing = GF(UnivariateRingZp64(7, "z")("1 + 3*z + z^2 + z^3"), "z")
// ring GF(7,3)[a,b,c]
implicit val ring = MultivariateRing(cfRing, Array("a", "b", "c"))

// parse poly in GF(7^3)[a,b,c] from string
// coefficients of polynomials in GF(7,3)[a,b,c] are elements
// of GF(7,3) that is polynomials over "z"
val poly = ring("1 - (1 - z^3) * a^6*b + (1 - 2*z) * c^33 + a^66")

//factorize poly
println(Factor(poly))
```



### Rational function arithmetic

Define a field of rational functions *Frac(Z[x,y,z])* and input some functions:

``` scala
// Frac(Z[x,y,z])
implicit val field = Frac(MultivariateRing(Z, Array("x", "y", "z")))

// parse some math expression from string
// it will be automatically reduced to a common denominator
// with the gcd being automatically cancelled
val expr1 = field("(x/y/(x - z) + (x + z)/(y - z))^2 - 1")

// do some math ops programmatically
val (x, y, z) = field("x", "y", "z")
val expr2 = expr1.pow(2) + x / y - z
```

Greatest common divisors of numerators and denominators are always cancelled automatically. 

Use ``Coder`` to parse more complicated expressions:

``` scala
// bind expr1 and expr2 to variables to use them further in parser
field.coder.bind("expr1", expr1)
field.coder.bind("expr2", expr2)

// parse some complicated expression from string
// it will be automatically reduced to a common denominator
// with the gcd being automatically cancelled
val expr3 = field(
  """
     expr1 / expr2 - (x*y - z)/(x-y)/expr1
     + x / expr2 - (x*z - y)/(x-y)/expr1/expr2
     + x^2*y^2 - z^3 * (x - y)^2
  """)

// export expression to string
println(field.stringify(expr3))

// take numerator and denominator
val num = expr3.numerator()
val den = expr3.denominator()
// common GCD is always cancelled automatically
assert( field.ring.gcd(num, den).isOne )
```

Compute unique factor decomposition of rational function:

``` scala
// compute unique factor decomposition of expression
val factors = field.factor(expr3)
println(field.stringify(factors))
```

### Ideals and Groebner bases

Construct some ideal and check its properties:

``` scala
// ring Z/17[x,y,z]
implicit val ring = MultivariateRingZp64(17, Array("x", "y", "z"))
val (x, y, z) = ring("x", "y", "z")

// create ideal with two generators using GREVLEX monomial order for underlying Groebner basis
val I = Ideal(ring, Seq(x.pow(2) + y.pow(12) - z, x.pow(2) * z + y.pow(2) - 1), GREVLEX)
// I is proper ideal
assert(I.isProper)

// get computed Groebner basis
val gb = I.groebnerBasis
println(gb)

// check some ideal properties
assert(I.dimension == 1)
assert(I.degree == 36)
```

Unions, intersections and quotients of ideals:

``` scala
// create another ideal with only one generator
val J = Ideal(ring, Seq(x.pow(4) * y.pow(4) + 1), GREVLEX)
// J is principal ideal
assert(J.isPrincipal)

val union = I union J
// union is zero dimensional ideal
assert(union.dimension == 0)

val intersection = I intersection J
// intersection is still 2-dimensional
assert(intersection.dimension == 2)

// yet another ideal
val K = Ideal(ring, Seq(z * x.pow(4) - z * y.pow(14) + y * z.pow(16), (x + y + z).pow(4)), GREVLEX)
// compute complicated quotient ideal
val quotient = (I * J * K) :/ times
assert(quotient == K) 
```


------------------------------------------------------------------------

Construct lexicographic Gröbner basis to solve a system of equations:

``` scala
// ring Q[a, b, c]
implicit val ring = MultivariateRing(Q, Array("x", "y", "z"))

// parse some polynomials from strings
val a = ring("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z")
val b = ring("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4")
val c = ring("8*x^3 + 12*y^3 + x*z^2 + 3")
val d = ring("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3")

// construct ideal with Groebner basis in LEX order
val ideal = Ideal(ring, Seq(a, b, c, d), LEX)
// it is very simple: <z^2, x, 1+4*y^3>
println(ideal)
```


### Programming

Implement generic function for solving linear Diophantine equations:

``` scala
/**
  * Solves equation \sum f_i s_i  = gcd(f_1, \dots, f_N) for given f_i and unknown s_i
  * @return a tuple (gcd, solution)
  */
def solveDiophantine[E](fi: Seq[E])(implicit ring: Ring[E]) =
  fi.foldLeft((ring(0), Seq.empty[E])) { case ((gcd, seq), f) =>
    val xgcd = ring.extendedGCD(gcd, f)
    (xgcd(0), seq.map(_ * xgcd(1)) :+ xgcd(2))
  }
```

Implement generic function for computing partial fraction decomposition:

``` scala
/** Computes partial fraction decomposition of given rational */
def apart[E](frac: Rational[E]) = {
  implicit val ring: Ring[E] = frac.ring
  val factors = ring.factor(frac.denominator).map {case (f, exp) => f.pow(exp)}
  val (gcd,  nums) = solveDiophantine(factors.map(frac.denominator / _))
  val (ints, rats) = (nums zip factors)
    .map { case (num, den) => Rational(frac.numerator * num, den * gcd) }
    .flatMap(_.normal)       // extract integral parts from fractions
    .partition(_.isIntegral) // separate integrals and fractions
  rats :+ ints.foldLeft(Rational(ring(0)))(_ + _)
}
```

Apply that function to elements of different rings:

``` scala
// partial fraction decomposition for rationals
// gives List(184/479, (-10)/13, 1/8, (-10)/47, 1)
val qFracs = apart( Q("1234213 / 2341352"))

// partial fraction decomposition for rational functions
val ufRing = Frac(UnivariateRingZp64(17, "x"))
// gives List(4/(16+x), 1/(10+x), 15/(1+x), (14*x)/(15+7*x+x^2))
val pFracs = apart( ufRing("1 / (3 - 3*x^2 - x^3 + x^5)") )
```

------------------------------------------------------------------------

Implement Lagrange method for univariate interpolation:

```latex
    p(x) = \sum_i p(x_i) \Pi_{j \ne i} \frac{x_{\phantom{i}} - x_j}{x_i -x_j}
```

``` scala
/** Lagrange polynomial interpolation formula */
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
```

Interpolate polynomial from *Frac(Z_{13}[a,b,c])[x]*:

``` scala
// coefficient ring Frac(Z/13[a,b,c])
val cfRing = Frac(MultivariateRingZp64(2, Array("a", "b", "c")))
val (a, b, c) = cfRing("a", "b", "c")

implicit val ring = UnivariateRing(cfRing, "x")
// interpolate with Lagrange formula
val data = Seq(a -> b, b -> c, c -> a)
val poly = interpolate(data)
assert(data.forall { case (p, v) => poly.eval(p) == v })
```


Highlighted benchmarks
----------------------

<img src="https://github.com/PoslavskySV/rings/blob/develop/doc/_static/gcd_nvars.png?raw=true" width="600">

Dependence of multivariate GCD performance on the number of variables. For details see [benchmarks](https://github.com/PoslavskySV/rings.benchmarks/tree/master/gcd)

<img src="https://github.com/PoslavskySV/rings/blob/develop/doc/_static/gcd_size.png?raw=true" width="600">

Dependence of multivariate GCD performance on the size of input polynomials. For details see [benchmarks](https://github.com/PoslavskySV/rings.benchmarks/tree/master/gcd)

<img src="https://github.com/PoslavskySV/rings/blob/develop/doc/_static/factor.png?raw=true" width="600">

Dependence of multivariate factorization performance on the number of variables. For details see [benchmarks](https://github.com/PoslavskySV/rings.benchmarks/tree/master/factor)

<img src="https://github.com/PoslavskySV/rings/blob/develop/doc/_static/factor_univar.png?raw=true" width="600">

Univariate factorization performance on polynomials of the form *(1 + \sum_{i = 1}^{i \leq deg} i \times x^i)* in *Z_{17}[x]*.


Index of algorithms implemented in Rings
----------------------------------------

The list of algorithms implemented in Rings (and references to the literature) can be found at [Rings RTD page](http://rings.readthedocs.io/en/latest/algorithms.html)


------------------------------------------------------------------------

License
-------

Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2>.
