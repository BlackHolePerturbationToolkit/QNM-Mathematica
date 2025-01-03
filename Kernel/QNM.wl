(* ::Package:: *)

(* ::Title:: *)
(*QNM*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["QNM`",
  {
  "SpinWeightedSpheroidalHarmonics`",
  "Teukolsky`",
  "Teukolsky`TeukolskyRadial`"
  }
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{QNMFrequency, QNMRadial, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, n, a] computes a quasinormal mode frequency.";


QNMRadial::usage = "QNMRadial[s, l, m, n, a, opts] computes the radial eigenfunction of a quasinormal mode.";


QNMRadialFunction::usage = "QNMRadialFunction[...] is an object representing a quasinormal mode solution to the radial Teukolsky equation.";


(* ::Subsection::Closed:: *)
(*Error messages*)


QNMFrequency::nointerp = "Interpolation data not available for s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::optx = "Unknown options in `1`";
QNMFrequency::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::findroot = "FindRoot failed to converge to the requested accuracy.";
QNMFrequency::cmplx = "Only real values of a are allowed, but a=`1` specified.";
QNMRadial::optx = "Unknown options in `1`";
QNMRadial::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMRadial::coords = "Coordinate options are either \"BL\", \"Boyer-Lindquist\", or \"Hyperboloidal, but got `1`";
QNMRadial::convergence = "Eigenvalue failed to converge to specified tolerance. Final value `1`";
QNMRadialFunction::dmval = "Radius `1` lies outside the computational domain.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Global Debug Flag*)


DEBUG=False;


(* ::Section::Closed:: *)
(*Hyperboloidal equations*)


(* ::Subsection::Closed:: *)
(*Operators*)


\[CapitalDelta][\[Rho]_,a_]:=1-2M \[Rho]+a^2 \[Rho]^2;
A[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=2 I \[Omega]-2(1+s)\[Rho]+2(I \[Omega](a^2-8M^2)+I m a +(s+3)M)\[Rho]^2+4(2 I \[Omega] M-1) a^2 \[Rho]^3;
B[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=(a^2-16M^2)\[Omega]^2+2(m a +2 I s M)\[Omega]+2(4(a^2-4M^2)M \[Omega]^2+(4 m a M-4I (s+2) M^2+I a^2)\[Omega]+ I m a+(s+1)M)\[Rho]+2(8 M^2 \[Omega]^2+6 I M \[Omega]-1)a^2 \[Rho]^2;
M\[Rho][\[Omega]_,s_Integer,m_Integer,a_]:=-\[Rho]^2\[CapitalDelta][\[Rho],a]R''[\[Rho]]+A[\[Omega],s,m ,a,\[Rho]]R'[\[Rho]]+B[\[Omega],s,m ,a,\[Rho]]R[\[Rho]];


(* ::Subsection::Closed:: *)
(*Horizon locations*)


M=1;


rp[a_, M_] := M+Sqrt[M^2-a^2];
rm[a_, M_] := M-Sqrt[M^2-a^2];


(* ::Subsection::Closed:: *)
(*Discretization*)


(* Define the grid *)
\[Rho]grid[a_, NN_Integer] := 1/rp[a, M] Reverse[1/2 (1+Cos[\[Pi] Subdivide[NN-1]])];


(* Get differentiation matrices based on the grid *)
Dgrid[\[Rho]grid_] := NDSolve`FiniteDifferenceDerivative[Derivative[1], \[Rho]grid, DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];
DDgrid[\[Rho]grid_] := NDSolve`FiniteDifferenceDerivative[Derivative[2], \[Rho]grid, DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];


(* ::Subsection::Closed:: *)
(*Radial Eigenvalues*)


(* Compute the eigenvalue of the Radial equation *)
\[Lambda]r[\[Omega]_, s_Integer, m_Integer, a_, NN_] :=
 Module[{\[Rho]grida, dd1, dd2, Mat, DiscretizationRules},
  (* Create grid and differentiation matrices *)
  \[Rho]grida = \[Rho]grid[a,NN];
  dd1 = Dgrid[\[Rho]grida];
  dd2 = DDgrid[\[Rho]grida];

  (* Define discretization in terms of grids *)
  DiscretizationRules={\[Rho]->\[Rho]grida, R''[\[Rho]]->dd2, R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[NN]};

  (* Evaluate discretized operator *)
  Mat = M\[Rho][\[Omega], s, m, a] /. DiscretizationRules;

  (* Return sorted eigenvalues *)
  SortBy[Eigenvalues[Mat], Norm]
]


(* ::Section::Closed:: *)
(*QNMFrequency*)


(* ::Subsection::Closed:: *)
(*Interpolation of tabulated QNM frequencies*)


$QNMInstallationDirectory = FileNameDrop[FindFile["QNM`"], -2];
$QNMDataDirectory = FileNameJoin[{$QNMInstallationDirectory, "Data"}];


ClearAll[QNMFrequencyInterpolation];

Options[QNMFrequencyInterpolation] = Options[Interpolation];

QNMFrequencyInterpolation[s_, l_, m_, n_, opts:OptionsPattern[]] := QNMFrequencyInterpolation[s, l, m, n] =
 Module[{h5file, dataset, data, ret},
  h5file = FileNameJoin[{$QNMDataDirectory, "QNM_s"<>ToString[s]<>".h5"}];
  dataset = "/l"<>ToString[l]<>"/m"<>ToString[m]<>"/n"<>ToString[n];
  Quiet[data = Import[h5file, {"Datasets", dataset}, "ComplexKeys"->{"r", "i"}];, {Import::general, Import::noelem, Import::nffil}];
  If[MatchQ[data, <|"a"->_, "omega"->_|>],
    ret = Interpolation[Transpose[Lookup[data, {"a","omega"}]], opts];,
    Message[QNMFrequency::nointerp, s, l, m, n];
    ret = $Failed;
  ];
  ret
];


(* ::Subsection::Closed:: *)
(*Calculate QNM frequency by finding zero of incidence amplitude*)


Options[QNMFrequencyInIncidenceAmplitude] = {"InitialGuess" -> Automatic};


QNMFrequencyInIncidenceAmplitude[s_Integer, l_Integer, m_Integer, n_, a_, OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, inInc},
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Check[QNMFrequencyInterpolation[s, l, m, n][a], 1, QNMFrequency::nointerp];
  ];
  inInc[\[Omega]_?NumericQ] := TeukolskyRadial[s, l, m, a, \[Omega], Method -> "MST"]["In"]["Amplitudes"]["Incidence"];
  \[Omega]QNM /. FindRoot[inInc[\[Omega]QNM], {\[Omega]QNM, \[Omega]guess}]
];


(* ::Subsection::Closed:: *)
(*Hyperboloidal method*)


(* ::Subsubsection::Closed:: *)
(*Eigenvalue difference*)


\[Delta]\[Lambda][\[Omega]guess_, s_Integer, l_Integer, m_Integer, a_, NN_Integer] :=
 Module[{\[Lambda]},
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]guess];
  -\[Lambda]r[\[Omega]guess, s, m, a, NN][[1]] - (\[Lambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2)
]


(* ::Subsubsection::Closed:: *)
(*Define function options*)


Options[QNMFrequencyHyperboloidal] = {"Tolerance"->10^-6, "NumPoints"->100, "InitialGuess"->Automatic};


(* ::Subsubsection::Closed:: *)
(*Newton Raphson solver*)


(*CP commment: the algorithm solves for the separation constant, while SpinWeightedSpheroidalEigenvalue gives the Angular eigenvalue: SeparationConstant= (SpinWeightedSpheroidalEigenvalue+2 m a \[Omega]guess- (a \[Omega]guess)^2) *)
QNMFrequencyHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, a_, opts:OptionsPattern[]] :=
 Module[
	{\[Lambda],\[Omega]guess,F,Fp,\[Epsilon],\[Gamma],\[Delta]\[Omega],tol,NN,MAXITS,count,monitor},
		(* Load option values *)
		tol = OptionValue["Tolerance"];
		NN = OptionValue["NumPoints"];
		(* Debug *)
		If[DEBUG,
		Print["Calculating QNMFrequency with tolerance ", N[tol], " for ", NN, " gridpoints."]
		];
		\[Omega]guess=OptionValue["InitialGuess"];
		If[\[Omega]guess === Automatic,
		  \[Omega]guess = SetPrecision[Check[QNMFrequencyInterpolation[s, l, m, n][a], 1, QNMFrequency::nointerp], Precision[a]];
		];
		If[DEBUG,
		monitor = PrintTemporary["Eigenvalue: ", Dynamic[\[Omega]guess]];
		];
		\[Delta]\[Omega]=1;
		\[Epsilon]=10^-8;(*better way?*)
		\[Gamma]=1;(*allow this to change?*)
		
		(* Iterate until tolerance is reached or maximum iterations *)
		MAXITS=1000;
		count = 0;
		Until[Norm[\[Delta]\[Omega]]<tol,
			count += 1;
			F=\[Delta]\[Lambda][\[Omega]guess, s,l, m, a,NN];
			Fp=(\[Delta]\[Lambda][\[Omega]guess+\[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-\[Epsilon], s,l, m, a,NN])/(2 \[Epsilon])+(\[Delta]\[Lambda][\[Omega]guess+I \[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-I \[Epsilon], s,l, m, a,NN])/(2 \[Epsilon] I);
			\[Omega]guess= \[Omega]guess-(\[Gamma] F)/Fp;
			\[Delta]\[Omega]= -(\[Gamma] F)/Fp;
						
			If[count > MAXITS, 
				Message[QNMRadial::convergence, \[Omega]guess];
				Return[$Failed];
			];
		];
		If[DEBUG,
		NotebookDelete[monitor];
		];
		
		(* Once result is obtained, return eigenvalue & separation constant *)
		\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
		<|"QNMfrequency"->\[Omega]guess,"SeparationConstant"->\[Lambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2|>
];


(* ::Subsection::Closed:: *)
(*QNMFrequency*)


SyntaxInformation[QNMFrequency] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[QNMFrequency] = {Method -> Automatic};


SetAttributes[QNMFrequency, {NumericFunction, Listable, NHoldAll}];


QNMFrequency[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] || !IntegerQ[n] || n < 0 :=
 (Message[QNMFrequency::params, s, l, m, n]; $Failed);


QNMFrequency[s_, l_, m_, n_, a_Complex, OptionsPattern[]] :=
 (Message[QNMFrequency::cmplx, a]; $Failed);


QNMFrequency[s_, l_, m_, n_, a_?InexactNumberQ, OptionsPattern[]] := 
 Module[{opts, \[Omega]ini, \[Omega]},
  Switch[OptionValue[Method],
    "Interpolation",
      \[Omega] = QNMFrequencyInterpolation[s, l, m, n][a],
    {"Interpolation", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyInterpolation]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyInterpolation[s, l, m, n, opts][a];,
    "IncidenceAmplitude",
      \[Omega] = QNMFrequencyInIncidenceAmplitude[s, l, m, n, a],
    {"IncidenceAmplitude", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyInIncidenceAmplitude]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyInIncidenceAmplitude[s, l, m, n, a, opts];,
    Automatic | "Hyperboloidal",
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a]["QNMfrequency"];,
    {"Hyperboloidal", Rule[_,_]...},
	  opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a, opts]["QNMfrequency"];,
    _,
      Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      \[Omega] = $Failed;
  ];
  \[Omega]
];


QNMFrequency /: N[QNMFrequency[s_, l_, m_, n_, a_?NumericQ, opts:OptionsPattern[]], Nopts___] :=
  QNMFrequency[s, l, m, n, N[a, Nopts], opts];


(* ::Section::Closed:: *)
(*QNMRadial*)


(* ::Subsection::Closed:: *)
(*QNMRadialHyperboloidal*)


Options[QNMRadialHyperboloidal] = {
  "Coordinates" -> "Hyperboloidal",
  "NumPoints" -> 100
};


QNMRadialHyperboloidal[s_, l_, m_, n_, a_, \[Omega]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], ev, ef, \[Rho]grida, dd1, dd2, DiscretizationRules, Mat, RadialFunction, h, h\[Phi], NN, coords},
  (* Load options values *)
  NN = OptionValue["NumPoints"];

  (* Check a valid Coordinates option has been specified *)
  coords = OptionValue["Coordinates"];
  If[!MemberQ[{"BL", "Boyer-Lindquist", "BoyerLindquist","Hyperboloidal"}, coords],
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Construct the discretized system *)
  \[Rho]grida = \[Rho]grid[a,NN];
  dd1 = Dgrid[\[Rho]grida];
  dd2 = DDgrid[\[Rho]grida];

  (* Define discretization in terms of grids *)
  DiscretizationRules = {\[Rho]->\[Rho]grida, R''[\[Rho]]->dd2, R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[NN]};

  (* Evaluate discretized operator *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  ev = -(\[Lambda] + 2 a m \[Omega] - (a \[Omega])^2);
  Mat = M\[Rho][\[Omega],s,m,a] - ev R[\[Rho]] /. DiscretizationRules;

  (* Calculate the eigenvectors *)
  ef = First[NullSpace[Mat]];

  Switch[coords,
  "Hyperboloidal",
    RadialFunction = Function[{r}, Evaluate[Interpolation[Transpose[{\[Rho]grida,ef/ef[[-1]]}]][1/r]]];,
  "BL" | "BoyerLindquist" | "Boyer\[Dash]Lindquist",
    RadialFunction = Function[{r}, Evaluate[
      With[{rp = rp[a, M], rm = rm[a, M]},
        h = (2 M rp )/(rp-rm) Log[r-rp]-(2 M rm )/(rp-rm) Log[r-rm]-r-4 M Log[r];
        h\[Phi] = a/(rp-rm) Log[(r-rp)/(r-rm)];
        Interpolation[Transpose[{\[Rho]grida,ef/ef[[-1]]}]][1/r] Exp[-I*\[Omega]*h+I*m*h\[Phi]]]]];,
  _,
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Return QNMRadialFunction *)
  QNMRadialFunction[<|"s" -> s, "l" -> l, "m" -> m, "n" -> n, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
    "Method" -> "SpectralHyperboloidal", "RadialFunction" -> RadialFunction,
    "Coordinates" -> coords, "Domain" -> {rp[a, 1], \[Infinity]}|>]
]


(* ::Subsection::Closed:: *)
(*QNMRadial*)


SyntaxInformation[QNMRadial] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[QNMRadial] = {
  Method -> Automatic,
  "Frequency" -> Automatic
};


SetAttributes[QNMRadial, {Listable, NHoldAll}];


QNMRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] || !IntegerQ[n] || n < 0 :=
 (Message[QNMRadial::params, s, l, m, n]; $Failed);


QNMRadial[s_, l_, m_, n_, a_Complex, OptionsPattern[]] :=
 (Message[QNMRadial::cmplx, a]; $Failed);


QNMRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a:(_?InexactNumberQ|0), OptionsPattern[]] :=
 Module[{\[Omega], opts, qnm},
  (* Get the frequency *)
  If[OptionValue["Frequency"] === Automatic,
    \[Omega] = QNMFrequency[s, l, m, n, a];,
    \[Omega] = OptionValue["Frequency"];
  ];

  Switch[OptionValue[Method],
    Automatic|"Hyperboloidal",
      qnm = QNMRadialHyperboloidal[s, l, m, n, a, \[Omega]],
    {"Hyperboloidal", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMRadialHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMRadial::optx, Method -> OptionValue[Method]];
      ];
      qnm = QNMRadialHyperboloidal[s, l, m, n, a, \[Omega], opts],
    _,
     Message[QNMRadial::optx, Method -> OptionValue[Method]];
     qnm = $Failed;
  ];
  qnm
]


(* ::Section::Closed:: *)
(*QNMRadialFunction*)


SyntaxInformation[QNMRadialFunction] =
 {"ArgumentsPattern" -> {_}};


SetAttributes[QNMRadialFunction, {NHoldAll}];


(* ::Subsection::Closed:: *)
(*Output format*)


QNMRadialFunction /:
MakeBoxes[qnmrf: QNMRadialFunction[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", assoc["s"]}], "  ",
                  BoxForm`SummaryItem[{"l: ", assoc["l"]}], "  ",
                  BoxForm`SummaryItem[{"m: ", assoc["m"]}], "  ",
                  BoxForm`SummaryItem[{"n: ", assoc["n"]}]}],
             BoxForm`SummaryItem[{"a: ", assoc["a"]}]};
  extended = {BoxForm`SummaryItem[{"Frequency: ", assoc["\[Omega]"]}],
              BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Coordinates: ", assoc["Coordinates"]}]
              };
  BoxForm`ArrangeSummaryBox[
    QNMRadialFunction,
    qnmrf,
    None,
    summary,
    extended,
    form
  ]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


QNMRadialFunction[assoc_][key_String] := assoc[key];


QNMRadialFunction[assoc_]["RadialFunction"] := Missing["KeyAbsent", "RadialFunction"];


Keys[m_QNMRadialFunction] ^:= DeleteCases[Join[Keys[m[[-1]]], {}], "RadialFunction"];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


outsideDomainQ[r_, rmin_, rmax_] := Min[r]<rmin || Max[r]>rmax;


QNMRadialFunction[assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[QNMRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[assoc["RadialFunction"][r], InterpolatingFunction::dmval]
 ];


Derivative[n_][QNMRadialFunction[assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[QNMRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[Derivative[n][assoc["RadialFunction"]][r], InterpolatingFunction::dmval]
 ];


(*Derivative[n_Integer/;n>1][qnmrf:(QNMRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_])][r0:(_?NumericQ|{_?NumericQ..})] :=
 Module[{Rderivs, R, r, i, res},
  Rderivs = D[R[r_], {r_, i_}] :> D[0(*FIXME*), {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. Rderivs, {R'[r], R[r]}, Simplify];, {i, 2, n}];
  res = Derivative[n][R][r] /. {
    R'[r] -> qnmrf'[r0],
    R[r] -> qnmrf[r0], r -> r0};
  Clear[Rderivs, i];
  Remove[R, r];
  res
];*)


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency, QNMRadial, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
