import sbt.Keys._

organization := "cc.redberry"
name := "rings-scaladsl"
version := "1.0-SNAPSHOT"
scalaVersion := "2.12.3"

resolvers += Resolver.mavenLocal

libraryDependencies += "cc.redberry" % "rings" % "1.0-SNAPSHOT"
libraryDependencies += "junit" % "junit" % "4.12" % Test
