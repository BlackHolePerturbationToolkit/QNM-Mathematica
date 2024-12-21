(* ::Package:: *)

(* ::Title:: *)
(*QNM*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection:: *)
(*BeginPackage*)


BeginPackage["QNM`",
  {
  "SpinWeightedSpheroidalHarmonics`",
  "Teukolsky`",
  "Teukolsky`TeukolskyRadial`"
  }
];


(* ::Subsection:: *)
(*Unprotect symbols*)


ClearAttributes[{QNMFrequency, QNMMode, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, n, a, opts] computes the fundamental quasinormal frequency following arXiv:2202.03837. Options are:
\t\[Bullet] \"Tolerance\" -> \!\(\*SuperscriptBox[\(10\), \(-6\)]\) (default): Tolerance for Newton Raphson solver
\t\[Bullet] \"Resolution\" -> 100 (default): Spectral resolution
Returns: <| \"QNMfrequency\" -> \[Omega], \"SeparationConstant\" -> \[ScriptCapitalS] |>";


QNMMode::usage = "QNMMode[s, l, m, n, a, opts] computes the fundamental quasinormal mode frequency \[Omega] and corresponding radial eigenfunction following arXiv:2202.03837. Options are:
\t\[Bullet] \"Tolerance\" -> \!\(\*SuperscriptBox[\(10\), \(-6\)]\) (default): Tolerance for Newton Raphson solver
\t\[Bullet] \"Resolution\" -> 100 (default): Spectral resolution
\t\[Bullet] \"Coordinates\" -> {\"BL\", \"Boyer-Lindquist\", \"Hyperboloidal (default) \"}: Choice of output coordinates
Returns: <| \"QNMfrequency\" -> \[Omega], \"Radial Function\" -> F |>";


QNMRadialFunction::usage = "QNMRadialFunction[...] is an object representing a quasinormal mode solution to the radial Teukolsky equation.";


(* ::Subsection:: *)
(*Error messages*)


QNMFrequency::nointerp = "Interpolation data not available for s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::real = "Only real values of a are allowed, but a=`1`.";
QNMMode::coords = "Coordinate options are either \"BL\", \"Boyer-Lindquist\", or \"Hyperboloidal, but got `1`";
QNMMode::convergence = "Eigenvalue failed to converge to specified tolerance. Final value `1`";
QNMRadialFunction::dmval = "Radius `1` lies outside the computational domain.";


(* ::Subsection:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Global Debug Flag*)


DEBUG=False;


(* ::Section::Closed:: *)
(*Interpolation of tabulated QNM frequencies*)


$QNMInstallationDirectory = FileNameDrop[FindFile["QNM`"], -2];
$QNMDataDirectory = FileNameJoin[{$QNMInstallationDirectory, "Data"}];


ClearAll[QNMFrequencyInterpolation];

Options[QNMFrequencyInterpolation] = Options[Interpolation];

QNMFrequencyInterpolation[s_, l_, m_, n_, opts:OptionsPattern[]] := QNMFrequencyInterpolation[s, l, m, n] =
 Module[{h5file, dataset, data, ret},
  h5file = FileNameJoin[{$QNMDataDirectory, "QNM_s"<>ToString[s]<>".h5"}];
  dataset = "/l"<>ToString[l]<>"/m"<>ToString[m]<>"/n"<>ToString[n];
  Quiet[data = Import[h5file, {"Datasets", dataset}, "ComplexKeys"->{"r", "i"}];, {Import::general, Import::noelem}];
  If[MatchQ[data, <|"a"->_, "omega"->_|>],
    ret = Interpolation[Transpose[Lookup[data, {"a","omega"}]], opts];,
    Message[QNMFrequency::nointerp, s, l, m, n];
    ret = $Failed;
  ];
  ret
];


(* ::Section::Closed:: *)
(*Calculate QNM frequency by finding zero of incidence amplitude*)


QNMFrequencyInIncidenceAmplitude[s_Integer, l_Integer, m_Integer, a_, \[Omega]guess_] :=
 Module[{\[Omega]QNM, inInc},
  inInc[\[Omega]_?NumericQ] := TeukolskyRadial[s, l, m, a, \[Omega], Method -> "MST"]["In"]["Amplitudes"]["Incidence"];
  \[Omega]QNM /. FindRoot[inInc[\[Omega]QNM], {\[Omega]QNM, \[Omega]guess}]
];


(* ::Section:: *)
(*Calculate QNM Frequency*)


(*CP comment: To do: currently just n=0. Generalise?*)


(* ::Subsection:: *)
(*Parameters*)


M=1;


(* ::Subsection:: *)
(*Operator Functions*)


\[CapitalDelta][\[Rho]_,a_]:=1-2M \[Rho]+a^2 \[Rho]^2;
A[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=2 I \[Omega]-2(1+s)\[Rho]+2(I \[Omega](a^2-8M^2)+I m a +(s+3)M)\[Rho]^2+4(2 I \[Omega] M-1) a^2 \[Rho]^3;
B[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=(a^2-16M^2)\[Omega]^2+2(m a +2 I s M)\[Omega]+2(4(a^2-4M^2)M \[Omega]^2+(4 m a M-4I (s+2) M^2+I a^2)\[Omega]+ I m a+(s+1)M)\[Rho]+2(8 M^2 \[Omega]^2+6 I M \[Omega]-1)a^2 \[Rho]^2;
M\[Rho][\[Omega]_,s_Integer,m_Integer,a_]:=-\[Rho]^2\[CapitalDelta][\[Rho],a]R''[\[Rho]]+A[\[Omega],s,m ,a,\[Rho]]R'[\[Rho]]+B[\[Omega],s,m ,a,\[Rho]]R[\[Rho]];


(* ::Subsection:: *)
(*Discretization*)


rmax[a_]:=If[a==0, 1/2,1/a^2 (1-Sqrt[1-a^2]) ];
(* Define the grid *)
\[Rho]grid[a_, NN_Integer]:=rmax[a] Sort[N[1/2 (1+Cos[\[Pi] Range[0,1,1/(NN-1)]])]];
(* Get differentiation matrices based on the grid *)
Dgrid[\[Rho]grid_]:= Dgrid[\[Rho]grid]= NDSolve`FiniteDifferenceDerivative[Derivative[1],{\[Rho]grid},
	DifferenceOrder->{"Pseudospectral"},PeriodicInterpolation->{False}]["DifferentiationMatrix"];
DDgrid[\[Rho]grid_]:= DDgrid[\[Rho]grid] = NDSolve`FiniteDifferenceDerivative[Derivative[2],{\[Rho]grid},
	DifferenceOrder->{"Pseudospectral"},PeriodicInterpolation->{False}]["DifferentiationMatrix"];


(* ::Subsection:: *)
(*Radial Eigenvalues*)


(* Compute the eigenvalue of the Radial equation *)
\[Lambda]r[\[Omega]guess_,s_Integer, m_Integer, a_, NN_]:=Module[
	{\[Rho]grida,dd1,dd2,Mat,DiscretizationRules},
		(* Call grids *)
		\[Rho]grida = \[Rho]grid[a,NN];
		dd1 = Dgrid[\[Rho]grida];
		dd2 = DDgrid[\[Rho]grida];
		(* Define discretization in terms of grids *)
		DiscretizationRules={\[Rho]->\[Rho]grida,R''[\[Rho]]->dd2, 
			R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[NN]};
		(* Evaluate discretized operator *)
		Mat=M\[Rho][\[Omega]guess,s,m,a]/.DiscretizationRules;
		(* Return sorted eigenvalues *)
		Eigenvalues[{Mat}]//SortBy[#,Norm]&
]


(* ::Subsection:: *)
(*Eigenvalue difference*)


(* CP comment: Computes the difference between the  radial and angular separation constant! *)
\[Delta]\[Lambda][\[Omega]guess_,s_Integer,l_Integer,m_Integer,a_,NN_Integer]:=Module[
	{\[CapitalLambda]},
		\[CapitalLambda]=SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
		SetPrecision[-\[Lambda]r[\[Omega]guess,s,m,a,NN][[1]]-(\[CapitalLambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2),50]
]


(* ::Subsection:: *)
(*Define function options*)


Options[QNMFrequency] = {"Tolerance"->10^-6, "Resolution"->100, "\[Omega]0"->1};
Default[QNMFrequency] = 0;


(* ::Section:: *)
(*Newton Raphson solver*)


(*CP commment: the algorithm solves for the separation constant, while SpinWeightedSpheroidalEigenvalue gives the Angular eigenvalue: SeparationConstant= (SpinWeightedSpheroidalEigenvalue+2 m a \[Omega]guess- (a \[Omega]guess)^2) *)
QNMFrequency[s_Integer,l_Integer,m_Integer,n_Integer,a_, opts:OptionsPattern[]]:=Module[
	{\[Lambda],\[Omega]guess,F,Fp,\[Epsilon],\[Gamma],\[Delta]\[Omega],tol,NN,MAXITS,count,monitor},
		(* Check for real spin *)
		If[!RealValuedNumberQ[a],
		Message[QNMFrequency::real, a];
		Return[$Failed];
		];
		
		(* Load option values *)
		tol = OptionValue["Tolerance"];
		NN = OptionValue["Resolution"];
		(* Debug *)
		If[DEBUG,
		Print["Calculating QNMFrequency with tolerance ", N[tol], " for ", NN, " gridpoints."]
		];
		\[Omega]guess=OptionValue["\[Omega]0"];
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
			F=SetPrecision[\[Delta]\[Lambda][\[Omega]guess, s,l, m, a,NN],50];
			Fp=SetPrecision[(\[Delta]\[Lambda][\[Omega]guess+\[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-\[Epsilon], s,l, m, a,NN])/(2 \[Epsilon])+(\[Delta]\[Lambda][\[Omega]guess+I \[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-I \[Epsilon], s,l, m, a,NN])/(2 \[Epsilon] I),50];
			\[Omega]guess= SetPrecision[\[Omega]guess-(\[Gamma] F)/Fp,50];
			\[Delta]\[Omega]= SetPrecision[-(\[Gamma] F)/Fp,50];
						
			If[count > MAXITS, 
				Message[QNMMode::convergence, \[Omega]guess];
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


(* ::Subsection:: *)
(*Overload for fewer arguments*)


QNMFrequency[s_Integer,l_Integer,m_Integer,a_,opts:OptionsPattern[]] := QNMFrequency[s,l,m,Default[QNMFrequency],a,opts];


(* ::Section:: *)
(*Calculate radial function*)


rp[a_,M_] := M+Sqrt[M^2-a^2];
rm[a_,M_] := M-Sqrt[M^2-a^2];


(* ::Subsection:: *)
(*Define function options*)


Options[QNMMode] = {"Tolerance"->10^-6, "Resolution"->100, "Coordinates"->"Hyperboloidal"};
Default[QNMMode] = 0;


(* ::Subsection:: *)
(*Use the eigenvalue to calculate the eigenfunction*)


(* Calculate the eigenvalue first, then solve for the radial profile *)
QNMMode[s_Integer,l_Integer,m_Integer,n_Integer,a_, opts:OptionsPattern[]]:=Module[
	{\[Omega],ev,ef,\[Rho]grida,dd1,dd2,DiscretizationRules,Mat,RadialFunction,Delta,
	DeltaTilde,h,h\[Phi],tol,NN,coords},
		(* Load options values *)
		tol = OptionValue["Tolerance"];
		NN = OptionValue["Resolution"];
		coords = OptionValue["Coordinates"];
		
		(* Debug *)
		If[DEBUG,
			Print["Computing QNMMode with tolerance ", N[tol], " resolution ", NN, " and coordinates ", coords];
		];
		
		If[coords != "BL" && coords != "Boyer-Lindquist" && coords != "BoyerLindquist" && coords != "Hyperboloidal",
		Message[QNMMode::coords, coords];
		Return[$Failed];
		];
		
		(* Get the frequency *)
		\[Omega]=QNMFrequency[s,l,m,0,a,"Tolerance"->tol, "Resolution"->NN]["QNMfrequency"];
		
		(* Call the discretized system *)
		\[Rho]grida = \[Rho]grid[a,NN];
		dd1 = Dgrid[\[Rho]grida];
		dd2 = DDgrid[\[Rho]grida];
		
		(* Define discretization in terms of grids *)
		DiscretizationRules={\[Rho]->\[Rho]grida,R''[\[Rho]]->dd2, 
			R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[NN]};
			
		(* Evaluate discretized operator *)
		Mat=M\[Rho][\[Omega],s,m,a]/.DiscretizationRules;
		
		(* Calculate the eigensystem given the eigenvalue *)
		{ev,ef}=Sort[Transpose[Eigensystem[{Mat}]], Norm[#1[[1]]]<Norm[#2[[1]]]&][[1]];
		
		Switch[coords,
		"Hyperboloidal",
			RadialFunction = Interpolation[Transpose[{\[Rho]grida,ef/ef[[-1]]}]];,
		"BL" | "BoyerLindquist" | "Boyer\[Dash]Lindquist" ,
			h=-(1/\[Rho]grida)+(2 M^2 ArcTan[(-M+1/\[Rho]grida)/Sqrt[a^2-M^2]])/Sqrt[a^2-M^2]+M Log[a^2+1/\[Rho]grida^2-(2 M)/\[Rho]grida]-4 M Log[1/\[Rho]grida];
			h\[Phi]=(a ArcTan[(-M+1/\[Rho]grida)/Sqrt[a^2-M^2]])/Sqrt[a^2-M^2];
			RadialFunction = Interpolation[Transpose[{\[Rho]grida,ef/ef[[-1]] Exp[-I*\[Omega]*h+I*m*h\[Phi]]}]];,
		_,
			Message[QNMMode::coords, coords];
			Return[$Failed];
		];

		(* Return association *)
		QNMRadialFunction[<|"s" -> s, "l" -> l, "m" -> m, "n" -> n, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> ev,
		  "Method" -> "SpectralHyperboloidal", "RadialFunction" -> RadialFunction,
		  "Coordinates" -> coords, "Domain" -> {rp[a, 1], \[Infinity]}
		  |>]
]


(* ::Subsection:: *)
(*Overload for fewer arguments*)


QNMMode[s_Integer,l_Integer,m_Integer, a_, opts:OptionsPattern[]] := QNMMode[s,l,m,Default[QNMMode],a,opts];


(* ::Section:: *)
(*QNMRadialFunction*)


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
              BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"],
              BoxForm`SummaryItem[{"Coordinates: ", assoc["Coordinates"]}]
              }]};
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


QNMRadialFunction[assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];


Keys[m_QNMRadialFunction] ^:= DeleteCases[Join[Keys[m[[-1]]], {}], "RadialFunction"];


(* ::Subsection:: *)
(*Numerical evaluation*)


SetAttributes[QNMRadialFunction, {NumericFunction}];


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


Derivative[n:1][QNMRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[QNMRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[Derivative[n][assoc["RadialFunction"]][r], InterpolatingFunction::dmval]
 ];


Derivative[n_Integer/;n>1][qnmrf:(QNMRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_])][r0:(_?NumericQ|{_?NumericQ..})] :=
 Module[{Rderivs, R, r, i, res},
  Rderivs = D[R[r_], {r_, i_}] :> D[0(*FIXME*), {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. Rderivs, {R'[r], R[r]}, Simplify];, {i, 2, n}];
  res = Derivative[n][R][r] /. {
    R'[r] -> qnmrf'[r0],
    R[r] -> qnmrf[r0], r -> r0};
  Clear[Rderivs, i];
  Remove[R, r];
  res
];


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency, QNMMode, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
