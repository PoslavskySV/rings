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

### Quick start


