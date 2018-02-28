.. role:: rings
  :class: rings

.. role:: feature
  :class: feature

.. |br| raw:: html

   <br/>


.. |Groebner| raw:: html

   Gr&ouml;bner



Rings: efficient Java/Scala library for polynomial rings
########################################################

|Rings| is an efficient implementation of univariate and multivariate polynomial algebra over arbitrary coefficient rings. It makes use of asymptotically fast algorithms for basic algebraic operations as well as for advanced methods like GCDs, factorization and |Groebner| bases. Performance achieved in |Rings| is comparable to well known solutions like Singular/NTL/FLINT/Maple/Mathematica.

The key features of |Rings| include:

 * :ref:`Polynomials → <ref-basics-polynomials>` |br| :feature:`Univariate and multivariate polynomials over arbitrary coefficient rings`

 * :ref:`Polynomial GCD → <ref-polynomial-methods>` |br| :feature:`Polynomial GCD over arbitrary coefficient domains`

 * :ref:`Univariate factorization → <ref-univariate-factorization>` |br| :feature:`Univariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
   
 * :ref:`Multivariate factorization → <ref-multivariate-factorization>` |br| :feature:`Multivariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
 
 * :ref:`Commutative algebra → <ref-rings>` |br| :feature:`Arbitrary rings, Galois fields, polynomial ideals etc`

 * :ref:`Ideals and Gröbner bases → <ref-ideals>` |br| :feature:`Polynomial ideals and efficient algorithms for Gröbner bases`

 * :ref:`Scala DSL → <ref-scala-dsl>` |br| :feature:`Powerful domain specific language in Scala`
   
 * :ref:`Fast → <ref-benchmarks>` |br| :feature:`Really fast library suitable for real-world computational challenges`


For a quick overview of what |Rings| can do proceed to :ref:`ref-quickstart` and try out |Rings|\ *.repl*.

Rings are hosted at GitHub: https://github.com/PoslavskySV/rings.

.. toctree::
   :caption: Documentation:
   :maxdepth: 3

   quickstart
   guide
   benchmarks
   algorithms



License
=======

Apache License, Version 2.0 http://www.apache.org/licenses/LICENSE-2.0.txt

