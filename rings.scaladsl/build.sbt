import sbt.Keys._

organization := "cc.redberry"
name := "rings-scaladsl"
version := "1.0-SNAPSHOT"
scalaVersion := "2.12.3"

resolvers += Resolver.mavenLocal

libraryDependencies ++= Seq(
  "cc.redberry" % "rings" % "1.0-SNAPSHOT",
  "junit" % "junit" % "4.12" % Test,
  "com.novocode" % "junit-interface" % "0.11" % Test exclude("junit", "junit-dep")
)