implicit val ring = MultivariateRing(Z, Array("x1", "x2", "x3", "x4", "x5", "x6", "x7"))

val poly1 = ring("1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7").pow(15) - 1
val poly2 = ring("1 + 3*x1*x2 + 5*x2*x3 + 7*x3*x4 + 9*x4*x5 + 11*x5*x6 + 13*x6*x7 + 15*x7*x1").pow(3) *
            ring("1 + 3*x1*x3 + 5*x2*x4 + 7*x3*x5 + 9*x6*x5 + 11*x7*x6 + 13*x6*x1 + 15*x7*x2").pow(3) *
            ring("1 + 3*x1*x4 + 5*x2*x5 + 7*x3*x6 + 9*x6*x7 + 11*x7*x1 + 13*x6*x2 + 15*x7*x3").pow(3) - 1
            

for (cfRing <- Seq(Z, Zp(2), Zp(524287))) {

   val tPoly1 = poly1.setRing(cfRing)
   val tPoly2 = poly2.setRing(cfRing)
   val tRing  = MultivariateRing(cfRing, ring.variables)   

   val (seconds1, factors1) = timing { tRing.factor(tPoly1) }
   println(s"Factor poly1 in ${cfRing}: ${seconds1}s, #factors: ${factors1.size}")

   val (seconds2, factors2) = timing { tRing.factor(tPoly2) }
   println(s"Factor poly2 in ${cfRing}: ${seconds2}s, #factors: ${factors2.size}")

}
