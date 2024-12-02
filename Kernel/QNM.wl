(* ::Package:: *)

(* ::Title:: *)
(*QNM*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection:: *)
(*BeginPackage*)


BeginPackage["QNM`",
  {
  "SpinWeightedSpheroidalHarmonics`"
  }
];


(* ::Subsection:: *)
(*Unprotect symbols*)


ClearAttributes[{QNMFrequency}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, n, a] computes the quasinormal mode frequncy";


(* ::Subsection:: *)
(*Error messages*)


QNMFrequency::real = "Only real values of a are allowed, but a=`1`.";


(* ::Subsection:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*Frequency*)


QNMFrequency[s_Integer, l_Integer, m_Integer, n_Integer, a_] :=
 Module[{\[Lambda], \[Omega]guess},
  If[!RealValuedNumberQ[a],
    Message[QNMFrequency::real, a];
    Return[$Failed];
  ];
  \[Omega]guess = 0.3 + 0.2 I;
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
  \[Lambda]
]


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
