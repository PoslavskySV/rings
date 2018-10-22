import sbt.Keys._

organization := "cc.redberry"

name := "rings.scaladsl"

version := "2.5.3"

scalaVersion := "2.12.3"

crossScalaVersions := Seq("2.11.11", "2.12.3")

moduleName := name.value

resolvers += Resolver.mavenLocal

libraryDependencies ++= Seq(
  "cc.redberry" % "rings" % "2.5.3",
  "junit" % "junit" % "4.12" % Test,
  "com.novocode" % "junit-interface" % "0.11" % Test exclude("junit", "junit-dep")
)

publishTo := Some(
  if (isSnapshot.value)
    Opts.resolver.sonatypeSnapshots
  else
    Opts.resolver.sonatypeStaging
)

import sbtrelease.ReleasePlugin.autoImport.ReleaseTransformations._

releaseCrossBuild := true // true if you cross-build the project for multiple Scala versions
releaseProcess := Seq[ReleaseStep](
  checkSnapshotDependencies,
  runClean,
  //runTest,
  releaseStepCommand("publishSigned")
)
