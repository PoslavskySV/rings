import java.io._
import scala.io._

val exp = 7
{// GCD in Z[X]

   implicit val ring = MultivariateRing(Z, Array("x1", "x2", "x3", "x4", "x5", "x6", "x7"))

   val a = ring("1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7").pow(exp) - 1
   val b = ring("1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7").pow(exp) + 1
   val g = ring("1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7").pow(exp) + 3

   val ag = a*g
   val bg = b*g

   val start1 = System.nanoTime()
   val r1 = ring.gcd(ag, bg)
   println(s"Done in Z: ${System.nanoTime() - start1}")
   assert (r1 % g == ring(0))

   val start2 = System.nanoTime()
   val r2 = ring.gcd(ag + 1, bg)
   println(s"Done in Z (trivial): ${System.nanoTime() - start2}")
}

{// GCD in Zp[X]

   implicit val ring = MultivariateRing(Zp(524287), Array("x1", "x2", "x3", "x4", "x5", "x6", "x7"))

   val a = ring("1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7").pow(exp) - 1
   val b = ring("1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7").pow(exp) + 1
   val g = ring("1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7").pow(exp) + 3

   val ag = a*g
   val bg = b*g

   val start1 = System.nanoTime()
   val r1 = ring.gcd(ag, bg)
   println(s"Done in Zp: ${System.nanoTime() - start1}")
   assert (r1 % g == ring(0))

   val start2 = System.nanoTime()
   val r2 = ring.gcd(ag + 1, bg)
   println(s"Done in Zp (trivial): ${System.nanoTime() - start2}")
}