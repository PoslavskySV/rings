[![Build Status](https://circleci.com/gh/PoslavskySV/r2.0x.svg?style=shield&circle-token=2db6894ad742a4f5147d1a3089ec91d4c637a7c3)](https://circleci.com/gh/PoslavskySV/r2.0x)

# Rings: efficient Java/Scala library for polynomial rings

`Rings` provides very efficient implementation of polynomial algebra for both univariate and multivariate cases over arbitrary coefficient rings. It uses asymptotically fast algorithms for basic algebraic operations as well as for advanced methods like GCDs and polynomial factorization and gives performance comparable or sometimes even better than well known solutions like `Singular`/`Maple`/`Mathematica`.

Highlighted features of `Rings` are:

 * univariate and multivariate polynomials over arbitrary coefficient rings
 * polynomial factorization (univariate and multivariate) and GCDs over arbitrary coefficient rings
 * efficient and arbitrary large Galois fields
 * really high performance


### Installation

`Rings` is available from Maven Central:

```xml
<dependency>
    <groupId>cc.redberry</groupId>
    <artifactId>rings</artifactId>
    <version>1.0</version>
</dependency>
```

For Scala interface:

```xml
<dependency>
    <groupId>cc.redberry</groupId>
    <artifactId>rings.scaladsl</artifactId>
    <version>1.0</version>
</dependency>
```

or via SBT:

```scala
libraryDependencies += "cc.redberry" % "rings.scaladsl" % "1.0"
```

### Examples


Polynomial factorization in Z/17[x]:

```scala
import cc.redberry.rings.poly.PolynomialMethods._
import cc.redberry.rings.scaladsl.Rings._
import cc.redberry.rings.scaladsl.implicits._

// Ring Z/17[x]
implicit val ring = UnivariateRingZp64(17, "x")

// parse univariate polynomials from string
val poly1 = ring("1 + 2*x + 3*x^2")
val poly2 = ring("4 + 5*x^5 + 6*x^6")
val poly3 = ring("7 + 8*x^9 + 9*x^10 + 10*x^11 + x^99")
val poly = poly1 * poly2 * poly3

// factorize polynomial
val factors = Factor(poly)
// there will be exactly 9 factors
assert(9 == factors.size())

println(s"factorization : ${ring show factors}")
```

One can use coefficient ring with arbitrary large characteristic. For example replacing `val ring = ...` in the above example with

```scala
// Z/1237940039285380274899124357
val coefficientRing = Zp(new BigInteger("1267650600228229401496703205653"))
// Z/1237940039285380274899124357[x]
implicit val ring = UnivariateRing(coefficientRing, "x")
```
will factor `poly` in Z/1267650600228229401496703205653[x] (the next prime to 2^100).
