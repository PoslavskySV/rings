========================================
Index of algorithms implemented in Rings
========================================



Univariate polynomials
======================

1. *Karatsuba multiplication* [vzGG03]_ (Sec. 8.1) used (in some adapted form) for multiplication of univariate polynomials: 

 - `UnivariatePolynomial.multiply <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariatePolynomial.java>`_
 - `UnivariatePolynomialZp64.multiply <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariatePolynomialZp64.java>`_
	 
2. *Half-GCD and Extended Half-GCD* [vzGG03]_ (Sec. 11) used (in some adapted form inspired by [NTL]_) for univariate GCD:

 - `UnivariateGCD.HalfGCD  <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_
 - `UnivariateGCD.ExtendedHalfGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_
 
3. *Subresultant polynomial remainder sequences* [GCL92]_ (Sec. 7.3):

 - `UnivariateGCD.SubresultantRemainders <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_

4. *Modular GCD in* :math:`Z[x]` *and* :math:`Q[x]` [vzGG03]_ (Sec. 6.7), small primes version is used:

 - `UnivariateGCD.ModularGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateGCD.java>`_

5. *Fast univariate division with Newton iteration* [vzGG03]_ (Sec. 9.1) used everywhere where multiple divisions (remainders) by the same divider are performed:

 - `UnivariateDivision.fastDivisionPreConditioning <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateDivision.java>`_
 - `UnivariateDivision.divideAndRemainderFast <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateDivision.java>`_
 
6. *Univariate square-free factorization in zero characteristic (Yun's algorithm)* [vzGG03]_ (Sec. 14.6):

 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
     
7. *Univariate square-free factorization in non-zero characteristic (Musser's algorithm)* [Mus71]_, [GCL92]_ (Sec. 8.3):

 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationMusser <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
 - `UnivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateSquareFreeFactorization.java>`_
 
8. *Distinct-degree factorization* [vzGG03]_ (Sec. 14.2) plain version and adapted version (which is actually used) with precomputed :math:`x`-powers:

 - `DistinctDegreeFactorization.DistinctDegreeFactorizationPlain <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_
 - `DistinctDegreeFactorization.DistinctDegreeFactorizationPrecomputedExponents <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_

9. *Shoup's baby-step giant-step algorithm for distinct-degree factorization* [Sho95]_ used for factorization over fields with large cardinality:

 - `DistinctDegreeFactorization.DistinctDegreeFactorizationShoup <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/DistinctDegreeFactorization.java>`_

10. *Univariate modular composition* plain algorithm with Horner schema:
 
 - `ModularComposition.compositionHorner <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/ModularComposition.java>`_

11. *Brent-Kung algorithm* for univariate modular composition [BK98]_, [Sho95]_:

 - `ModularComposition.compositionBrentKung <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/ModularComposition.java>`_

12. *Cantor-Zassenhaus algorithm (equal-degree splitting)* [vzGG03]_ (Sec. 14.3) both for odd and even characteristic:

 - `EqualDegreeFactorization.CantorZassenhaus <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/EqualDegreeFactorization.java>`_

13. *Univaraite linear p-adic Hensel lifting* [GCL92]_ (Sec. 6.5):

 - `univar.HenselLifting.createLinearLift <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_
 - `univar.HenselLifting.liftFactorization <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_

14. *Univaraite quadratic p-adic Hensel lifting* [vzGG03]_ (Sec. 15.4-15.5):

 - `univar.HenselLifting.createQuadraticLift <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_
 - `univar.HenselLifting.liftFactorization <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/HenselLifting.java>`_

15. *Univariate polynomial factorization over finite fields* Musser's square free factorization -> distinct-degree factorization (either x-powers or Shoup's algorithm) -> Cantor-Zassenhaus equal-degree factorization:

 - `UnivariateFactorization.FactorInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_

16. *Univariate polynomial factorization over Z and Q* factorization modulo small prime -> Hensel lifting (adaptive linear/quadratic):

 - `UnivariateFactorization.FactorInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_
 - `UnivariateFactorization.FactorInQ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateFactorization.java>`_

17. *Univariate irreducibility test* [vzGG03]_ (Sec. 14.9):

 - `IrreduciblePolynomials.irreducibleQ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/IrreduciblePolynomials.java>`_

18. *Ben-Or’s generation of irreducible polynomials* [vzGG03]_ (Sec. 14.9):

 - `IrreduciblePolynomials.randomIrreduciblePolynomial <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/IrreduciblePolynomials.java>`_

19. *Univariate polynomial interpolation* Lagrange and Newton methods:

 - `UnivariateInterpolation <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/univar/UnivariateInterpolation.java>`_


Multivariate polynomials
========================


20. *Brown GCD over finite fields* [Bro71]_, [GCL92]_ (Sec. 7.4), [Yan09]_:

 - `MultivariateGCD.BrownGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

21. *Zippel's sparse GCD over finite fields* [Zip79]_, [Zip93]_, [dKMW05]_, [Yan09]_ both for monic (with fast Vandermonde systems) and non-monic (LINZIP) cases:

 - `MultivariateGCD.ZippelGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

22. *Extended Zassenhaus GCD (EZ-GCD) over finite fields* [GCL92]_ (Sec. 7.6), [MY73]_:

 - `MultivariateGCD.EZGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

23. *Enhanced Extended Zassenhaus GCD (EEZ-GCD) over finite fields* [Wan80]_:

 - `MultivariateGCD.EEZGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

24. *Modular GCD over Z with sparse interpolation* [Zip79]_, [Zip93]_, [dKMW05]_ (the same interpolation techniques as in `ZippelGCD` is used):

 - `MultivariateGCD.ModularGCD <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

25. *Kaltofen's & Monagan's generic modular GCD* [KM99]_ used for computing multivariate GCD over finite fields of very small cardinality

 - `MultivariateGCD.ModularGCDInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateGCD.java>`_

26. *Multivariate square-free factorization in zero characteristic (Yun's algorithm)* [Lee13]_:

 - `MultivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_

27. *Multivariate square-free factorization in non-zero characteristic (Musser's algorithm)* [Mus71]_, [GCL92]_ (Sec. 8.3):

- `MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_
- `MultivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateSquareFreeFactorization.java>`_

28. *Bernardin's fast dense multivariate Hensel lifting* [Ber99]_, [Lee13]_ both for bivariate case (original Bernardin's paper) and multivariate case (Lee thesis) and both with and without precomputed leading coefficients:

- `multivar.HenselLifting <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/HenselLifting.java>`_

29. *Fast dense bivariate factorization with recombination* [Ber99]_, [Lee13]_:

- `MultivariateFactorization.bivariateDenseFactorSquareFreeInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_
- `MultivariateFactorization.bivariateDenseFactorSquareFreeInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

30. *Kaltofen's multivariate factorization in finite fields* [Kal85]_, [Lee13]_; modified version of original Kaltofen's algorithm for leading coefficient precomputation with square-free decomposition (instead of distinct variables decomposition) due to Lee is used; further adaptations are made to work in finite fields of very small cardinality; the resulting algorithm is close to [Lee13]_, but at the same time has many differences in details:

- `MultivariateFactorization.factorInGF <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

31. *Kaltofen's multivariate factorization Z* [Kal85]_, [Lee13]_ (the same modifications as for finite field algorithm are made):

- `MultivariateFactorization.factorInZ <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateFactorization.java>`_

32. *Multivariate polynomial interpolation with Newton method*:

- `MultivariateInterpolation <https://github.com/PoslavskySV/rings/tree/develop/rings/src/main/java/cc/redberry/rings/poly/multivar/MultivariateInterpolation.java>`_
 

References
==========

.. [vzGG03] J. von zur Gathen and J. Gerhard. Modern computer algebra (2. ed.). Cambridge University Press, 2003.

.. [NTL] V. Shoup. NTL: A library for doing number theory. www.shoup.net/ntl

.. [GCL92] K. O. Geddes, S. R. Czapor, G. Labahn. Algorithms for Computer Algebra. 1992.

.. [Mus71] D.R. Musser,  Algorithms for polynomial factorization, Ph.D. Thesis, University of Wisconsin, 1971.

.. [Sho95] V. Shoup. A new polynomial factorization algorithm and its implementation. J. Symb. Comput., 20(4):363–397, 1995.

.. [BK98] R.P. Brent and H.T. Kung. Fast algorithms for manipulating formal power series. J. Assoc. Comput. Math. 25:581-595, 1978

.. [Bro71] W. S. Brown. On Euclid’s algorithm and the computation of polynomial greatest common divisors. J. ACM, 18(4):478–504, 1971.

.. [Zip79] R. E. Zippel. Probabilistic algorithms for sparse polynomials. In Proceedings of the International Symposiumon on Symbolic and Algebraic Computation, EUROSAM ’79, pages 216–226, London, UK, UK, 1979. Springer-Verlag.

.. [Zip93] R. E. Zippel. Effective Polynomial Computation. Kluwer International Series in Engineering and Computer Science. Kluwer Academic Publishers, 1993.

.. [dKMW05] J. de Kleine, M. B. Monagan, A. D. Wittkopf. Algorithms for the Non-monic Case of the Sparse Modular GCD Algorithm. Proceeding of ISSAC ’05, ACM Press, pp. 124-131 , 2005.

.. [Yan09] S. Yang. Computing the greatest common divisor of multivariate polynomials over finite fields. Master’s thesis, Simon Fraser University, 2009.

.. [MY73] J. Moses and D.Y.Y.Yun, "The EZGCD Algorithm," pp. 159-166 in Proc. ACM Annual Conference, (1973).

.. [Wan80] P.S. Wang, "The EEZ-GCD Algorithm," ACM SIGSAMBull., 14 pp. 50-60 (1980).

.. [KM99] E. Kaltofen, M. B. Monagan. On the Genericity of the Modular Polynomial GCD Algorithm. Proceeding of ISSAC ’99, ACM Press, 59-66, 1999.

.. [Ber99] L. Bernardin. Factorization of Multivariate Polynomials over Finite Fields. PhD thesis, ETH Zu ̈rich, 1999.

.. [Lee13] M. M.-D. Lee, Factorization of multivariate polynomials,  Ph.D. thesis, University of Kaiserslautern, 2013

.. [Kal85] E. Kaltofen. Sparse Hensel lifting. In EUROCAL 85 European Conf. Comput. Algebra Proc. Vol. 2, pages 4–17, 1985.


