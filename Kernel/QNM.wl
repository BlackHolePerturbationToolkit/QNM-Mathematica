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


QNMRadial::usage = "QNMRadial[s, l, m, n, a] computes the radial eigenfunction of a quasinormal mode.";


QNMRadialFunction::usage = "QNMRadialFunction[...] is an object representing a quasinormal mode solution to the radial Teukolsky equation.";


(* ::Subsection::Closed:: *)
(*Error messages*)


QNMFrequency::nointerp = "Interpolation data not available for s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::optx = "Unknown options in `1`.";
QNMFrequency::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::findroot = "FindRoot failed to converge to the requested accuracy.";
QNMFrequency::cmplx = "Only real values of a are allowed, but a=`1` specified.";
QNMFrequency::nokerr = "Method \"`1`\" only supported for Schwarzschild spacetime, but a=`2` specified.";
QNMFrequency::acc = "Accuracy of the calculated quasinormal mode frequency `1` is lower than that of the initial guess `2`.";
QNMFrequency::noacc = "Failed estimating accuracy of the calculated quasinormal mode frequency `1` compared to the initial guess `2`.";
QNMRadial::optx = "Unknown options in `1`.";
QNMRadial::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMRadial::coords = "Coordinate options are either \"BL\", \"Boyer-Lindquist\", or \"Hyperboloidal\", but got `1`.";
QNMRadial::convergence = "Eigenvalue failed to converge to specified tolerance. Final value `1`.";
QNMRadialFunction::dmval = "Radius `1` lies outside the computational domain.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Global Debug Flag*)


DEBUG=False;


(* ::Section::Closed:: *)
(*Utilities*)


(* ::Subsection::Closed:: *)
(*Horizon locations*)


M=1;


rp[a_, M_] := M+Sqrt[M^2-a^2];
rm[a_, M_] := M-Sqrt[M^2-a^2];


(* ::Subsection::Closed:: *)
(*Discretised Teukolsky operator on a hyperboloidal slice*)


\[ScriptCapitalM][s_, m_, a_, \[Omega]_, \[Lambda]_, n_] :=
 Module[{\[Rho], \[CapitalDelta], R, dR, d2R, A, B, M=1},
  \[Rho] = 1/rp[a, M] Reverse[1/2 (1+Cos[\[Pi] Subdivide[n-1]])];
  \[CapitalDelta] = 1-2M \[Rho]+a^2 \[Rho]^2;
  R = IdentityMatrix[n];
  dR  = NDSolve`FiniteDifferenceDerivative[Derivative[1], \[Rho], DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];
  d2R = NDSolve`FiniteDifferenceDerivative[Derivative[2], \[Rho], DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];
  A = 2 I \[Omega]-2(1+s)\[Rho]+2(I \[Omega](a^2-8M^2)+I m a +(s+3)M)\[Rho]^2+4(2 I \[Omega] M-1) a^2 \[Rho]^3;
  B = (a^2-16M^2)\[Omega]^2+2(m a +2 I s M)\[Omega]+2(4(a^2-4M^2)M \[Omega]^2+(4 m a M-4I (s+2) M^2+I a^2)\[Omega]+ I m a+(s+1)M)\[Rho]+2(8 M^2 \[Omega]^2+6 I M \[Omega]-1)a^2 \[Rho]^2;
  -\[Rho]^2 \[CapitalDelta] d2R + A dR + (B+(\[Lambda] + 2 a m \[Omega] - (a \[Omega])^2)) R
]


(* ::Subsection::Closed:: *)
(*Chebyshev interpolation*)


(* Convert from values at Chebyshev nodes to Chebyshev coefficients *)
chebCoeffs[values_] :=
 Module[{coeffs},
  coeffs = Sqrt[2/(Length[values]-1)] FourierDCT[values, 1];
  coeffs[[{1,-1}]] /= 2;
  coeffs
];


(* Construct a Chebyshev InterpolatingFunction *)
(* Useful references on the structure of an InterpolatingFunction:
   https://mathematica.stackexchange.com/a/28341
   https://mathematica.stackexchange.com/a/98349
   https://mathematica.stackexchange.com/a/114065 *)
chebInterp[data_, domain_] :=
 Module[{order, coeffs},
  order = Length[data] - 1;
  coeffs = chebCoeffs[data];
  InterpolatingFunction[
    {domain}    (* Interpolation domain *),
    {5          (* InterpolatingFunction version *),
     1,         (* Bit field: always 1 for Chebyshev *)
     order,     (* Max derivative order in data *)
     {2},       (* Domain grid size *)
     {order+1}, (* Interpolation order + 1 *)
     0,         (* Not a derivative of an existing InterpolatingFunction *)
     0,         (* Non-periodic interpolation *)
     0, 0,
     Automatic, (* Extrapolation handler *)
     {}, {}, False},
    {domain},   (* Domain grid *)
    {data, coeffs},
    {{{{1}, {1,2}, {1,2}}, {Automatic, ChebyshevT, ChebyshevT}}}
  ]
];


(* ::Subsection::Closed:: *)
(*Functions for Leaver's method*)


\[Beta][s_] := s^2 -1;
k1[m_, s_] := 1/2 Abs[m-s];
k2[m_, s_] := 1/2 Abs[m+s];


(* ::Subsubsection::Closed:: *)
(*Schwarzschild*)


(* Functions for the Continued Fraction used in Leaver's method for Schwarzschild *)
\[Alpha][i_, \[Omega]_]  := i^2 - 4 M I \[Omega] i + 2 i - 4 M I \[Omega] + 1;
\[Delta][i_, \[Omega]_,s_, l_] := -2 i^2 - (2 - 16 M I \[Omega])i + 32 M^2 \[Omega]^2 + 8 M I \[Omega] - l(l + 1) + \[Beta][s];
\[Gamma][i_, \[Omega]_, s_] := i^2 - 8 M I \[Omega] i - 16 M^2 \[Omega]^2 - \[Beta][s] - 1;


(* ::Subsubsection::Closed:: *)
(*Kerr*)


b[a_] := Sqrt[4 M^2 - 4 a^2];


(* Functions for the Continued Fraction used in Leaver's method for Kerr *)
\[Alpha]freq[i_,\[Omega]_, s_, m_, a_]:= i^2 + (2-s-2 M I \[Omega] -2 I/b[a] (2 M^2 \[Omega] - a m))i + 1 - s - 2 M I \[Omega] -2 I/b[a] (2 M^2 \[Omega] - a m);
\[Beta]freq[i_,\[Omega]_, Alm_, s_, m_, a_]:= -2 i^2 +(-2 + 2(4 M + b[a]) I \[Omega] + 4 I/b[a] (2 M^2 \[Omega] - a m)) i + (16 M^2 + 4 M b[a] - a^2) \[Omega]^2 - s - 1 - 2 a m \[Omega] - Alm +(4 M + b[a]) I \[Omega] + (8 M \[Omega] + 2 I)/b[a] (2 M^2 \[Omega] - a m);
\[Gamma]freq[i_, \[Omega]_, s_, m_, a_]:= i^2 + (s - 6 M I \[Omega] - 2 I/b[a] (2 M^2 \[Omega] - a m)) i - 8 M^2 \[Omega]^2 - 4 M I \[Omega] s - 8 M \[Omega]/b[a] (2 M^2 \[Omega] - a m);

\[Alpha]ang[i_, \[Omega]_, s_, m_, a_]:= -2 (i + 1) (i + 2 k1[m, s] + 1);
\[Beta]ang[i_, \[Omega]_, Alm_, s_, m_, a_]:= i (i - 1) + 2 i (k1[m, s] + k2[m, s] + 1 - 2 a \[Omega]) - (2 a \[Omega] (2 k1[m, s] + s + 1) - (k1[m, s] + k2[m, s]) (k1[m, s] + k2[m, s] + 1)) - (a^2 \[Omega]^2 + s(s+1) + Alm);
\[Gamma]ang[i_, \[Omega]_, s_, m_, a_]:= 2 a \[Omega] (i + k1[m, s] + k2[m, s] + s);


(* ::Section::Closed:: *)
(*QNMFrequency*)


prec[a_] :=
 Which[
   a != 0,
   Precision[a],
   a == 0 && Precision[a] == MachinePrecision,
   MachinePrecision,
   a == 0,
   Accuracy[a]
 ];


(* ::Subsection::Closed:: *)
(*Asymptotic approximations*)


(* ::Text:: *)
(*Multiple asymptotic expansions are used for Schwarzschild, with different expansions providing better approximations for different values of l and n.*)


(* ::Subsubsection::Closed:: *)
(*Schwarzschild asymptotic expansion for l >> n, l>> 1*)


(* ::Text:: *)
(*This expansion is taken from Dolan, Ottewill, Classical and Quantum Gravity, Vol. 26, 2009.*)


Schwarzfinit1[s_, l_, n_]:= ((I*(n+1/2)*(\[Beta][s]^2/27 + (\[Beta][s]*(1100*(n+1/2)^2 - 2719))/46656 + (11273136*(n+1/2)^4 - 52753800*(n+1/2)^2 + 66480535)/2902376448))/(l+1/2)^4 + 
   (-(\[Beta][s]^2/27) + (\[Beta][s]*(204*(n+1/2)^2 + 211))/3888 + (854160*(n+1/2)^4 - 1664760*(n+1/2)^2 - 776939)/40310784)/(l+1/2)^3 - (I*(n+1/2)*(\[Beta][s]/9 + (235*(n+1/2)^2)/3888 - 1415/15552))/(l+1/2)^2 + 
   (\[Beta][s]/3 - (5*(n+1/2)^2)/36 - 115/432)/(l+1/2) + (l+1/2) - I*(n+1/2))/(Sqrt[27]*M);


(* ::Subsubsection::Closed:: *)
(*Schwarzschild low order asymptotic expansion for n>>l, n>>1*)


(* ::Text:: *)
(*The expansion is taken from Casals, Dolan, Ottewill, Wardell, Phys. Rev. D, Vol. 88, 2013*)


Schwarzfinit2[n_] := Log[3]/(8 \[Pi] M) - I (n+1/2)/(4 M);


(* ::Subsubsection::Closed:: *)
(*Kerr*)


(* ::Text:: *)
(*For Kerr, an initial guess is needed for both the frequency and the spheroidal eigenvalue Subscript[A, lm].*)


Kerrfinit[s_, l_, m_, n_, a_] :=
 Module[{b, \[CapitalDelta], \[Mu], Eikonal, Rp, \[CapitalOmega]r, \[CapitalOmega]i, finit},
  \[CapitalDelta][r_] := r^2 -2 M r + a^2;
  \[Mu] = m/(l+1/2);(* Useful parameter *)
  Eikonal[rp_]:= 2(rp/M)^4(rp/M - 3)^2 + 4 (rp/M)^2((1 - \[Mu]^2)(rp/M)^2 - 2(rp/M) - 3(1 - \[Mu]^2))(a/M)^2 + (1 - \[Mu]^2) ( (2 - \[Mu]^2) (rp/M)^2 + 2 (2 + \[Mu]^2) (rp/M) + (2 - \[Mu]^2)) (a/M)^4;

  Rp = FindRoot[Eikonal[rp]==0, {rp, 3}] [[1]][[2]]; (* FIXME: This always produces a machine-precision result. Maybe use Solve instead? *)

  \[CapitalOmega]r = -If[\[Mu] == 0, \[Pi]/2 Sqrt[\[CapitalDelta][Rp]]/((Rp^2 + a^2) EllipticE[(a^2 \[CapitalDelta][Rp])/((Rp^2 + a^2)^2)] ) , (M - Rp) \[Mu] a / ((Rp - 3M)Rp^2 + (Rp + M) a^2)];
  \[CapitalOmega]i = \[CapitalDelta][Rp](Sqrt[4(6Rp^2 \[CapitalOmega]r^2 -1) + 2a^2\[CapitalOmega]r^2(3 - \[Mu]^2)])/(2 Rp^4 \[CapitalOmega]r - 4 a M Rp \[Mu] + a^2 Rp \[CapitalOmega]r(Rp(3 - \[Mu]^2) + 2 M(1+\[Mu]^2)) + a^4\[CapitalOmega]r(1-\[Mu]^2));

  finit = (l + 1/2) Abs[\[CapitalOmega]r] - I (n + 1/2) Abs[\[CapitalOmega]i]
];


KerrAinit[s_, l_, m_, n_, a_] := (l+1/2)^2 - (a Kerrfinit[s, l, m, n, a])^2 /2 (1 - (m/(l+1/2))^2);


(* ::Subsection::Closed:: *)
(*Interpolation of tabulated QNM frequencies*)


$QNMInstallationDirectory = FileNameDrop[FindFile["QNM`"], -2];
$QNMDataDirectory = FileNameJoin[{$QNMInstallationDirectory, "Data"}];


ClearAll[QNMFrequencyInterpolation];

Options[QNMFrequencyInterpolation] = Options[Interpolation];

QNMData[file_, dataset_] := QNMData[file, dataset] = Quiet[Import[file, {"Datasets", dataset}, "ComplexKeys"->{"r", "i"}];, {Import::general, Import::noelem, Import::nffil}];

QNMFrequencyInterpolation[s_, l_, m_, n_, opts:OptionsPattern[]] :=
 Module[{h5file, dataset, data, ret},
  h5file = FileNameJoin[{$QNMDataDirectory, "QNM_s"<>ToString[s]<>".h5"}];
  dataset = "/l"<>ToString[l]<>"/m"<>ToString[m]<>"/n"<>ToString[n];
  data = QNMData[h5file, dataset];
  If[MatchQ[data, <|"a"->_, "omega"->_|>],
    ret = Interpolation[Transpose[Lookup[data, {"a","omega"}]], opts];,
    Message[QNMFrequency::nointerp, s, l, m, n];
    ret = Function[{a}, $Failed];
  ];
  ret
];


(* ::Subsection::Closed:: *)
(*Calculate QNM frequency by finding zero of incidence amplitude*)


Options[QNMFrequencyInIncidenceAmplitude] = {"InitialGuess" -> Automatic};


QNMFrequencyInIncidenceAmplitude[s_Integer, l_Integer, m_Integer, n_, a_, OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, inInc, prec = prec[a]},
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Quiet[QNMFrequencyInterpolation[s, l, m, n][a], QNMFrequency::nointerp];
  ];
  If[\[Omega]guess == $Failed, 
    Message[QNMFrequency::nointerp, s, l, m, n];
    If[a==0,
      Which[
        l <= Max[n,2], 
        \[Omega]guess = SpectralInitialGuess[s, l, n];,
        l>n,
        \[Omega]guess = Schwarzfinit1[s, l, n];
      ];,
      \[Omega]guess = Kerrfinit[s, l, m, n, a];
    ];
  ];
  inInc[\[Omega]_?NumericQ] := inInc[\[Omega]] = TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega], Method -> "MST"]["In"]["Amplitudes"]["Incidence"];
  \[Omega]QNM /. Quiet[Check[FindRoot[inInc[\[Omega]QNM], {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum], FindRoot::nlnum]
];


(* ::Subsection::Closed:: *)
(*Ripley Hyperboloidal method*)


Options[QNMFrequencyHyperboloidal] = {"NumPoints" -> 32, "InitialGuess" -> Automatic, "AccuracyCheck" -> True};


QNMFrequencyHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, a_, opts:OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, \[Omega], \[Delta]\[Lambda], numpoints, prec = prec[a]},
  numpoints = OptionValue["NumPoints"];
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Quiet[QNMFrequencyInterpolation[s, l, m, n][a], QNMFrequency::nointerp];
  ];
  If[\[Omega]guess == $Failed, 
    Message[QNMFrequency::nointerp, s, l, m, n];
    If[a==0,
      Which[
        l <= Max[n,2], 
        \[Omega]guess = SpectralInitialGuess[s, l, n];,
        l>n,
        \[Omega]guess = Schwarzfinit1[s, l, n];
      ];,
      \[Omega]guess = Kerrfinit[s, l, m, n, a];
    ];
  ];
  \[Delta]\[Lambda][\[Omega]_?NumericQ] := \[Delta]\[Lambda][\[Omega]] = Module[{\[Lambda], Mat},
    \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
    Mat = \[ScriptCapitalM][s, m, a, \[Omega], \[Lambda], numpoints];
    Eigenvalues[Mat, -1]
  ];
  \[Omega]QNM = \[Omega] /. FindRoot[\[Delta]\[Lambda][\[Omega]], {\[Omega], SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec];

  If[OptionValue["AccuracyCheck"] == True &&
     Quiet[Check[Abs[TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega]QNM, Method -> "MST"]["In"]["Amplitudes"]["Incidence"] / 
               TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega]guess, Method -> "MST"]["In"]["Amplitudes"]["Incidence"]] > 10^3,
           Message[QNMFrequency::noacc, \[Omega]QNM, \[Omega]guess]; False], All, QNMFrequency::noacc],
    Message[QNMFrequency::acc, \[Omega]QNM, \[Omega]guess];
  ];

  \[Omega]QNM
];


(* ::Subsection::Closed:: *)
(*Ansorg-Macedo hyperboloidal method*)


(* ::Text:: *)
(*Spectral method using hyperboloidal slicing. This technique was adapted from a Reissner-Nordstrom version kindly provided by Rodrigo Macedo, Queen Mary University of London. The technique is outlined in detail in Ansorg, Macedo, Phys. Rev. D, Vol. 93, 2016*)


SpectralInitialGuess[s_, l_, n_] :=
 Module[{L1a, L1b, L2a, L2b, L2c, Prec, Ndiv, ndiv, \[Sigma], \[Sigma]0, \[Sigma]1, \[CapitalDelta]\[Sigma], x, Id, Z, c, d\[Sigma], D\[Sigma], D2\[Sigma], L1, L2, M, Eigens, Filtered, sPattern},
  (* Operators appearing in the wave equation in coordinates (\[Tau], \[Sigma]), excluding radial derivatives *)
  L1a = -((2 \[Sigma])/(\[Sigma]+1));
  L1b = (1-2\[Sigma]^2)/(\[Sigma]+1);
  L2a = (-l(l+1) - \[Sigma](1-s^2))/(\[Sigma]+1);
  L2b = (\[Sigma](2-3\[Sigma]))/(\[Sigma]+1);
  L2c = (\[Sigma]^2 (1-\[Sigma]))/(\[Sigma]+1);

  Prec = 100; (*Numerical precision, machine precision found to be insufficient*)
	  		(* Make optional for users to set? Requires more testing as well *)
  If[n > 2, Ndiv = 6n, Ndiv = 20]; (* Subdivision of radial domain \[Sigma] \[Element] [0,1] used for discretization of radial derivs *)
  ndiv = Ndiv + 1; (* Matrix dimension *)

  \[Sigma]0 = 0;
  \[Sigma]1 = 1;
  \[CapitalDelta]\[Sigma] = \[Sigma]1-\[Sigma]0;

  (* The method for discretizing the radial derivatives is detailed in Spectral Methods in Matlab, Lloyd N. Trefethen *)
  x[i_]:= Cos[(i \[Pi])/Ndiv]; (* Chebyshev points *)
  \[Sigma] = N[Table[\[Sigma]0 + \[CapitalDelta]\[Sigma]  (1+x[i])/2, {i, 0, Ndiv}], Prec]; (* Radial coordinate grid *)

  Id = IdentityMatrix[ndiv];
  Z = ConstantArray[0, {ndiv, ndiv}]; (* Zero matrix *)
  c[i_] := If[i(i-Ndiv) ==0, 2, 1];

  (* Radial derivative using Chebyshev-Lobatto method *)
  d\[Sigma][i_, j_] := 2/\[CapitalDelta]\[Sigma] Which[i==j && i == 0, (2 Ndiv^2+1)/6, i==j && i == Ndiv, -((2 Ndiv^2+1)/6), i==j,  -x[i]/(2(1-x[i]^2)) , i!=j, c[i]/c[j] (-1)^(i+j)/(x[i]-x[j])] ;

  (* First deriv matrix *)
  D\[Sigma]=N[Table[d\[Sigma][i, j], {i, 0, Ndiv}, {j, 0, Ndiv}], Prec];
  (* Second deriv matrix *)
  D2\[Sigma] = D\[Sigma] . D\[Sigma];

  (* Full opertors in wave eq *)
  L1 = L1b D\[Sigma] + L1a Id;
  L2 = L2c D2\[Sigma] + L2b D\[Sigma] + L2a Id;

  M = ArrayFlatten[{{Z, Id}, {L2, L1}}];

  (* Solve for eigenvalues *)
  Eigens = Eigenvalues[M] ;

  sPattern = Which[s==0,x_/;(Im[x]==0. ), 
    Abs[s]==1, x_/;(Im[x]==0. ) ,
    Abs[s]==2,  x_/;(Im[x]==0. && Abs[Re[x]+4] >0.2 )]; 

  (* Remove all values with Im[\[Omega]] = 0, which actually correspond to the purely imaginary roots (see multiplication by I/4 below).
     Also remove those with Im[\[Omega]] > 0, as these are the corresponding QNMs in the 3rd quadrant. *)
  Filtered = Reverse[DeleteCases[Eigens,x_/;(Im[x]>=0. )]];

  (*Pick out eigenvalue corresponding to the desired overtone *)
  (*Would be nice to save this and use the other elements for initial seeds for other values of n *)
  (*However, the accuracy decreases as one goes down the list *)
  Which[Abs[s]==2 && l==2 && n>8, I/4 Filtered[[n]], Abs[s]==2 && l==2 && n==8, (0.4615178773933189 10^-15 - I 3.9999999999996)/2, True, I/4 Filtered[[n+1]]]
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
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a];,
    {"Hyperboloidal", Rule[_,_]...},
	  opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a, opts];,
    "Ansorg-Macedo",
      If[a!=0,
        Message[QNMFrequency::nokerr, "Ansorg-Macedo", a];
        \[Omega] = $Failed;
        ,
        \[Omega] = SpectralInitialGuess[s, l, n];
      ];,
    "Large-l Asymptotic",
      If[a == 0,
        \[Omega] = Schwarzfinit1[s, l, n];,
        \[Omega] = Kerrfinit[s, l, m, n, a];
      ];,
    "Large-n Asymptotic",
      If[a!=0,
        Message[QNMFrequency::nokerr, "Large-n Asymptotic", a];
        \[Omega] = $Failed;
        ,
        \[Omega] = Schwarzfinit2[n];
      ];,
    _,
      Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      \[Omega] = $Failed;
  ];
  \[Omega]
];


QNMFrequency[s_, l_, m_, n_, 0, opts___] := QNMFrequency[s, l, m, n, 0.0, opts];


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
 Module[{\[Lambda], ef, ns, Mat, RadialFunction, h, h\[Phi], numpoints, coords},
  (* Load options values *)
  numpoints = OptionValue["NumPoints"];

  (* Check a valid Coordinates option has been specified *)
  coords = OptionValue["Coordinates"];
  If[!MemberQ[{"BL", "Boyer-Lindquist", "BoyerLindquist","Hyperboloidal"}, coords],
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Construct the discretized system *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  Mat = \[ScriptCapitalM][s, m, a, \[Omega], \[Lambda], numpoints];

  (* Calculate the eigenvectors *)
  ns = NullSpace[Mat];
  If[Length[ns]==0,
    ef = Eigensystem[Mat][[2, -1]];,
    ef = First[ns];
  ];

  Switch[coords,
  "Hyperboloidal",
    RadialFunction = Function[{r}, Evaluate[chebInterp[Reverse[ef/ef[[-1]]], {0, 1/rp[a, M]}][1/r]]];,
  "BL" | "BoyerLindquist" | "Boyer\[Dash]Lindquist",
    RadialFunction = Function[{r}, Evaluate[
      With[{rp = rp[a, M], rm = rm[a, M]},
        h = (2 M rp )/(rp-rm) Log[r-rp]-(2 M rm )/(rp-rm) Log[r-rm]-r-4 M Log[r];
        h\[Phi] = a/(rp-rm) Log[(r-rp)/(r-rm)];
        chebInterp[Reverse[ef/ef[[-1]]], {0, 1/rp[a, M]}][1/r] Exp[-I*\[Omega]*h+I*m*h\[Phi]]]]];,
  _,
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Return QNMRadialFunction *)
  QNMRadialFunction[<|"s" -> s, "l" -> l, "m" -> m, "n" -> n, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
    "Method" -> "SpectralHyperboloidal", "Amplitudes" -> <|"\[ScriptCapitalH]" -> ef[[-1]]/ef[[-1]], "\[ScriptCapitalI]" -> ef[[1]]/ef[[-1]]|>, "RadialFunction" -> RadialFunction,
    "Coordinates" -> coords, "Domain" -> {rp[a, M], \[Infinity]}|>]
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


QNMRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_?InexactNumberQ, OptionsPattern[]] :=
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


QNMRadial[s_, l_, m_, n_, 0, opts___] := QNMRadial[s, l, m, n, 0.0, opts];


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
