import sbt.Keys._

organization := "cc.redberry"
name := "rings.scaladsl"
version := "2.0"
scalaVersion := "2.12.3"
moduleName := name.value

resolvers += Resolver.mavenLocal

libraryDependencies ++= Seq(
  "cc.redberry" % "rings" % "2.0",
  "junit" % "junit" % "4.12" % Test,
  "com.novocode" % "junit-interface" % "0.11" % Test exclude("junit", "junit-dep")
)

crossScalaVersions := Seq("2.11.11", "2.12.3")