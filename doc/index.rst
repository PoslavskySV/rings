.. role:: rings
  :class: rings

.. role:: feature
  :class: feature

.. |br| raw:: html

   <br/>



Rings: efficient Java/Scala library for polynomial rings
########################################################

|Rings| is an efficient implementation of univariate and multivariate polynomial algebra over arbitrary coefficient rings. It makes use of asymptotically fast algorithms for basic algebraic operations as well as for advanced methods like GCDs and polynomial factorization. Performance achieved in |Rings| is comparable to well known solutions like Singular/NTL/FLINT/Maple/Mathematica.

The key features of |Rings| include:

 * :ref:`Polynomials → <ref-basics-polynomials>` |br| :feature:`Univariate and multivariate polynomials over arbitrary coefficient rings`

 * :ref:`Polynomial GCD → <ref-polynomial-methods>` |br| :feature:`Polynomial GCD over arbitrary coefficient domains`

 * :ref:`Univariate factorization → <ref-univariate-factorization>` |br| :feature:`Univariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
   
 * :ref:`Multivariate factorization → <ref-multivariate-factorization>` |br| :feature:`Multivariate polynomial factorization over arbitrary finite fields,` :math:`Z` :feature:`and` :math:`Q`
 
 * :ref:`Algebra → <ref-rings>` |br| :feature:`Arbitrary rings, Galois fields etc`

 * :ref:`Scala DSL → <ref-scala-dsl>` |br| :feature:`Powerful domain specific language in Scala`
   
 * :ref:`Fast → <ref-some-benchamrks>` |br| :feature:`Really fast library suitable for real-world computational challenges`


For a quick overview of what |Rings| can do proceed to :ref:`ref-quickstart` and try out |Rings|\ *.repl*.

.. toctree::
   :caption: Documentation:
   :maxdepth: 3

   quickstart
   guide
   algorithms


License
=======

Apache License, Version 2.0 http://www.apache.org/licenses/LICENSE-2.0.txt

