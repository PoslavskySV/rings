package cc.redberry.rings.scaladsl

import org.junit.Test

/**
  *
  */
class Examples3 {

  @Test
  def testNumberField1: Unit = {
    implicit val rationals = Q
    val minimalPolynomial = UnivariatePolynomial(3, 0, 0, 1)
    implicit val numberField = AlgebraicNumberField(minimalPolynomial, "a")

    println(numberField.minimalPolynomial(numberField("1 + a + a^2")))

    assert(numberField.degree() == 3)
    implicit val ring = UnivariateRing(numberField, "x")

    val poly = ring("(x^2 - a*x + a) * (x + a)")
    val factors = ring factor poly
    println(ring stringify factors)
  }

  @Test
  def testNumberField2: Unit = {
    implicit val rationals = Q
    val mPoly1 = UnivariatePolynomial(3, 0, 0, 1)

    implicit var cfField = MultipleFieldExtension(AlgebraicNumberField(mPoly1, "r1"))

    val mPoly2 = UnivariatePolynomial(2, 3, 4, 1)(cfField)
    cfField = cfField.joinAlgebraicElement(mPoly2, "r2")

    val mPoly3 = UnivariatePolynomial(cfField(2), cfField("r1 + r2"), cfField(4), cfField(1))(cfField)
    cfField = cfField.joinAlgebraicElement(mPoly3, "r3")

    println(cfField)
    println(cfField.degree())

    // very fast
    val simpleCfField = cfField.getSimpleExtension("A")
    simpleCfField.coder.bind("r1", cfField.getGeneratorRep(0))
    simpleCfField.coder.bind("r2", cfField.getGeneratorRep(1))
    simpleCfField.coder.bind("r3", cfField.getGeneratorRep(2))

    val psRing = MultivariateRing(simpleCfField, Array("x", "y", "z"))
    val sPoly = psRing("((x - r1 - r2) * (y - r1 - r3) * (z - r2 - r3))^12 - 1")
    println(sPoly.degree())

    //  very slow
//    val pmRing = MultivariateRing(cfField, Array("x", "y", "z"))
//    val mPoly = pmRing("((x - r1 - r2) * (y - r1 - r3) * (z - r2 - r3)) - 1")
//    println(mPoly.degree())
  }

  @Test
  def testComplexNumbers: Unit = {
    implicit val ring = MultivariateRing(GaussianRationals, Array("x", "y", "z"))
    val poly = ring("((x^2 - i)*(y^2 + i)*(z^2 - 2)*(i*x + i*y + i*z - 1))^3 - 1")
    val factors = ring.factor(poly)
    println(ring stringify factors)
    assert(2 == factors.size())
  }
}
