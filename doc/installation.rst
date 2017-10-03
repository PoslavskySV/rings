.. _ref-installation:

============
Installation
============


**Rings** are written in Java with DSL written in Scala. The Java library is available under ``cc.redberry.rings`` artifact from Maven Central:

.. code-block:: xml

	<dependency>
	    <groupId>cc.redberry</groupId>
	    <artifactId>rings</artifactId>
	    <version>1.0</version>
	</dependency>


Scala interface with DSL is available with both Maven:

.. code-block:: xml

	<dependency>
	    <groupId>cc.redberry</groupId>
	    <artifactId>rings.scaladsl</artifactId>
	    <version>1.0</version>
	</dependency>


and SBT:

.. code-block:: scala

	libraryDependencies += "cc.redberry" % "rings.scaladsl" % "1.0"


Sources are avalable at GitHub: https://github.com/PoslavskySV/rings