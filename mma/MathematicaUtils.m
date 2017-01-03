(****** utility functions *******)

IsSquareFree[poly_, modulus_] := FreeQ[PolynomialGCD[poly, D[poly, x], Modulus -> modulus], x]

$LC[poly_] := Coefficient[poly, x^(Exponent[poly, x])]

MonicPoly[poly_, modulus_] := Expand[ PowerMod[ $LC[poly], -1, modulus] poly, Modulus -> modulus]

RandomPoly[degree_, bound_] := Plus @@ Table[RandomInteger[{-bound, bound}] x^i, {i, 0, degree}];


(****** decode/endode factorization *******)

$RemoveSpaces[any_String] := StringReplace[any, {" " -> ""}];
$Encode[expr_] := ToString[expr, InputForm] // $RemoveSpaces

$EncodeSingleFactor[factor_] := Block[{f, exp},
  If[
    Head[factor] === Power,
    exp = factor[[2]]; f = factor[[1]];,
    exp = 1; f = factor;
  ];
  Return[$Encode[Join[{exp}, CoefficientList[f, x]]]]
];

$DecodeSingleFactor[factorStr_String] := Block[{list, f, exp},
  list = ToExpression /@ StringSplit[factorStr, ","];
  Return[(Plus @@ Table[list[[i]] x^(i - 2), {i, 2, Length[list]}]) ^ list[[1]] ];
];

EncodeFactorization[data_ModFactorization] := EncodeModFactorization @@ data;
EncodeModFactorization[poly_, factors_, modulus_] := Block[{exp, f, list},
  If[
    Head[poly] === Power,
    exp = poly[[2]]; f = poly[[1]];,
    exp = 1; f = poly
  ];
  list = Flatten@{$Encode[modulus], $EncodeSingleFactor @ poly, $EncodeSingleFactor /@ factors};
  StringReplace[StringRiffle[list, "|"] , {"{" -> "", "}" -> ""}] // $RemoveSpaces
];

DecodeModFactorization[factorizationStr_String] := Block[{list, modulus, poly, factors},
  list = StringSplit[factorizationStr, "|"];
  modulus = ToExpression[list[[1]]];
  list = $DecodeSingleFactor /@ list[[2 ;;]];
  Return[ModFactorization[list[[1]], list[[2 ;;]], modulus]];
];

(****** distinct-degree factorization *******)

$DistinctDegreeFactors[factorization_, modulus_] := Block[{exponents, dFactors, dFactor, fct},
  If[
    Head[factorization] === Times,
    fct = List @@ factorization,
    fct = {factorization}
  ];

  exponents = Exponent[fct, x] // DeleteDuplicates;

  dFactors = {};
  Do[
    AppendTo[dFactors, Expand[Times @@ Select[fct, Exponent[#, x] == exp &], Modulus -> modulus]];
    , {exp, exponents}];
  Return[dFactors];
];

FactorDistinctDegree[poly_, modulus_] := Block[{sqFactors, factors},
  sqFactors = FactorSquareFree[poly, Modulus -> modulus];
  If[
    Head[sqFactors] === Times,
    sqFactors = List @@ sqFactors,
    sqFactors = {sqFactors}
  ];
  factors = Flatten@(
    If[
      Head[#] === Power,
      $DistinctDegreeFactors[Factor[#[[1]], Modulus -> modulus], modulus]^#[[2]],
      $DistinctDegreeFactors[Factor[#, Modulus -> modulus], modulus]
    ]& /@ sqFactors
  );
  Return[factors];
];

RandomDistinctDegreeFactorization[maxDegree_, minNBase_, maxNBase_, bound_, primeBound_, seed_ : RandomInteger[2^32]] := Block[{nBase, poly, modulus, factors, sqFactors},
  SeedRandom[seed];
  While[True,
    nBase = RandomInteger[{minNBase, maxNBase}];
    poly = Times @@ Table[RandomPoly[RandomInteger[{1, maxDegree}], bound], {k, 1, nBase}];
    modulus = Prime[RandomInteger[{1, primeBound}]];
    poly = Expand[poly, Modulus -> modulus];
    If[
      poly == 0 || Head[poly] === Integer,
      Continue[]
    ];
    poly = MonicPoly[poly, modulus];
    Return[ModFactorization[poly, FactorDistinctDegree[poly, modulus], modulus]];
  ]
];


GenerateRandomDistinctFactorizationData[fileName_, nEntries_, initialSeed_, maxDegree_, minNBase_, maxNBase_, bound_, primeBound_ ] := Block[{out, tmp, timings},
  PrintTemporary[Dynamic[i]];
  SeedRandom[initialSeed];
  out = OpenWrite[fileName <> ".txt", BinaryFormat -> True, FormatType -> OutputForm, PageWidth -> Infinity];
  timings = {};
  Do[
    tmp = RandomDistinctDegreeFactorization[maxDegree, minNBase, maxNBase, bound, primeBound];
    AppendTo[timings, Timing[Factor[tmp[[1]], Modulus->tmp[[-1]] ]][[1]] ];
    tmp = EncodeFactorization[tmp];
    Write[out, tmp];
    , {i, nEntries}
  ];
  Close[out];

  CreateArchive[fileName <> ".txt", fileName <> ".gz"];
  DeleteFile[fileName <> ".txt"];
  Return[timings]
];

On[Assert];


(* Test decode/encode *)
Block[{i, tmp},
  For[i = 0, i < 100, ++i,
    tmp = RandomDistinctDegreeFactorization[4, 3, 6, 1000, 100];
    Assert[tmp[[1]] - Expand[Times @@ tmp[[2]], Modulus -> tmp[[-1]]] == 0];
    Assert[(tmp - DecodeModFactorization[EncodeFactorization[tmp]]) == 0];
  ];
]