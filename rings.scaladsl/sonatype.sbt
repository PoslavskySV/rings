
// Your profile name of the sonatype account. The default is the same with the organization value
sonatypeProfileName := "cc.redberry"

// To sync with Maven central, you need to supply the following information:
publishMavenStyle := true

// License of your choice
licenses := Seq("APL2" -> url("http://www.apache.org/licenses/LICENSE-2.0.txt"))
homepage := Some(url("https://github.com/PoslavskySV/rings"))
scmInfo := Some(
  ScmInfo(
    url("https://github.com/PoslavskySV/rings"),
    "scm:git@github.com:PoslavskySV/rings.git"
  )
)
developers := List(
  Developer(
    id = "PoslavskySV",
    name = "Stanislav Poslavsky",
    email = "stvlpos@mail.ru",
    url = url("http://redberry.cc")
  )
)
