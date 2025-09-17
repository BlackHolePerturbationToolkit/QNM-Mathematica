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
QNMFrequency::allln = "Method `1` can only be used with l=All and n=All to compute a set of frequencies but l=`2` and n=`3` specified.";
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
(*Discretised 2D Teukolsky operator on a hyperboloidal slice*)


\[ScriptCapitalM]2D[s_, m_, a_, nx_, nz_] :=
 Module[{\[Chi], z0, z1, \[CapitalDelta]z, z, x, Id, Zero, Dx, D2x, D\[Sigma], D2\[Sigma], L1\[Sigma]\[Sigma], L1\[Sigma], L10, L1xx, L1x, L2\[Sigma], L20, W, L1, L2,
         rh, \[Kappa], \[Stigma], \[Lambda], \[Delta]1, \[Delta]2, \[Sigma], X},
  (* Constants *)
  rh = rp[a, M];
  \[Kappa] = a/rh;
  \[Stigma] = s;
  \[Lambda] = M;
  \[Delta]1 = Abs[m-s];
  \[Delta]2 = Abs[m+s];

  (* Radial and angular coordinates *)
  z0 = 0;
  z1 = 1;
  \[CapitalDelta]z = z1-z0;
  \[Chi] = N@Reverse[Cos[\[Pi] Subdivide[nz-1]]];
  z = z0+1/2 \[CapitalDelta]z*(1+\[Chi]);
  x = N@Reverse[Cos[\[Pi] Subdivide[nx-1]]];

  (* Discretised operators *)
  Id = IdentityMatrix[nx*nz];
  Zero = ConstantArray[0, {nx*nz, nx*nz}];
  Dx  = NDSolve`FiniteDifferenceDerivative[Derivative[0,1], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D2x = NDSolve`FiniteDifferenceDerivative[Derivative[0,2], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D\[Sigma]  = NDSolve`FiniteDifferenceDerivative[Derivative[1,0], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D2\[Sigma] = NDSolve`FiniteDifferenceDerivative[Derivative[2,0], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];

  (* Coordinates on the full (flattened) 2D grid *)
  {\[Sigma], X} = {Flatten[Outer[#1&, z, x]], Flatten[Outer[#2&, z, x]]};

  L1xx = 1-X^2;
  L1x = \[Delta]1-\[Delta]2-X (2+\[Delta]1+\[Delta]2);
  
  L1\[Sigma]\[Sigma] = (-1+\[Sigma]) \[Sigma]^2 (-1+\[Kappa]^2 \[Sigma]);
  L1\[Sigma] = \[Sigma] (4 \[Kappa]^2 \[Sigma]^2+2 (1+\[Stigma])-\[Sigma] (2 I m \[Kappa]+(1+\[Kappa]^2) (3+\[Stigma])));
  L10 = 1/2 (-m^2-\[Delta]2-\[Delta]1 (1+\[Delta]2)-4 I m \[Kappa] \[Sigma]+4 \[Kappa]^2 \[Sigma]^2-2 (1+\[Kappa]^2) \[Sigma] (1+\[Stigma])+\[Stigma] (2+\[Stigma]));
  
  L2\[Sigma] = (2 rh (1+\[Sigma]^2 (-2+2 \[Kappa]^4 (-1+\[Sigma])+\[Kappa]^2 (-3+2 \[Sigma]))))/\[Lambda];
  L20 = (2 rh (3 (\[Kappa]^2+\[Kappa]^4) \[Sigma]^2-I m (\[Kappa]+2 \[Kappa] (1+\[Kappa]^2) \[Sigma])+\[Stigma]+\[Kappa] (-I X+\[Kappa]) \[Stigma]-\[Sigma] (2+\[Stigma]+\[Kappa]^4 (2+\[Stigma])+\[Kappa]^2 (3+2 \[Stigma]))))/\[Lambda];
  
  W = (rh^2 (4 (1+\[Sigma])+\[Kappa]^2 (7+X^2+4 \[Kappa]^2+4 (2+2 \[Kappa]^2+\[Kappa]^4) \[Sigma]-4 (1+\[Kappa]^2)^2 \[Sigma]^2)))/\[Lambda]^2;
  
  L1 = L1\[Sigma]\[Sigma] D2\[Sigma] + L1\[Sigma] D\[Sigma] + L10 Id + L1xx D2x + L1x Dx;
  L2 = L2\[Sigma] D\[Sigma] + L20 Id;
  ArrayFlatten[{{Zero, Id}, {L1/W, L2/W}}]
];


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


(* ::Subsection::Closed:: *)
(*Continued fraction*)


(* ::Text:: *)
(*This is preferred over Mathematica's ContinedFractionK function as we can get an error estimate on the result using this function.*)


CF[a_, b_, {n_, n0_}] := 
  Module[{A, B, ak, bk, res = Indeterminate, j = n0},
   A[n0 - 2] = 1;
   B[n0 - 2] = 0;
   ak[k_] := ak[k] = (a /. n -> k);
   bk[k_] := bk[k] = (b /. n -> k);
   A[n0 - 1] = 0(*bk[n0-1]*);
   B[n0 - 1] = 1;
   A[k_] := A[k] = bk[k] A[k - 1] + ak[k] A[k - 2];
   B[k_] := B[k] = bk[k] B[k - 1] + ak[k] B[k - 2];
   While[Quiet[Check[Abs[1-(A[j-1]/B[j-1])/(A[j]/B[j])], j--; False, General::munfl], General::munfl] > 10^(1-Precision[A[j]/B[j]]), j++];
   res = A[j]/B[j];
   Clear[A, B, ak, bk];
   res
];


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

QNMData[file_, dataset_] := QNMData[file, dataset] = Quiet[Import[file, {"Datasets", dataset}, "ComplexKeys"->{"r", "i"}], {Import::general, Import::noelem, Import::nffil}];

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
(*Assaad-Macedo hyperboloidal method in Kerr*)


(* ::Text:: *)
(*Spectral method using hyperboloidal slicing in Kerr. This technique was adapted from a version kindly provided by Rodrigo Macedo. The technique is outlined in detail in Assaad and Macedo, arXiv:2506.04326.*)


Options[SpectralInitialGuessKerr] = {"NumAngularPoints" -> 10, "NumRadialPoints" -> 10};


SpectralInitialGuessKerr[s_, m_, a_, opts:OptionsPattern[]] :=
 Module[{nAng, nRad},
  nAng = OptionValue["NumAngularPoints"];
  nRad = OptionValue["NumRadialPoints"];
  
  ReverseSortBy[I Eigenvalues[\[ScriptCapitalM]2D[s, m, a, nAng, nRad]], Im]
];


(* ::Subsection::Closed:: *)
(*Leaver method*)


(* ::Text:: *)
(*Calculation based on the method in Leaver, Proc. R. Soc. Lond. A. 402, 1985, Nollert, Phys. Rev. D, Vol. 47, 1993 as well as the code provided online by Emanuele Berti, https://pages.jh.edu/~eberti2/ringdown/.*)


(* Equations of which the QNMs and spheroidal eigenvalues are roots *)
Leaver[\[Omega]_?NumericQ, s_?IntegerQ, l_?IntegerQ, nInv_?IntegerQ] :=
 Module[{n},
  \[Delta][nInv, \[Omega], s, l] + ContinuedFractionK[-\[Alpha][nInv-n, \[Omega]] \[Gamma][nInv-n+1, \[Omega], s], \[Delta][nInv-n,\[Omega],s, l], {n,1,nInv}] + CF[-\[Alpha][n-1, \[Omega]] \[Gamma][n, \[Omega], s], \[Delta][n, \[Omega], s, l], {n, nInv+1}]
];
  
Leaver31[\[Omega]_?NumericQ, s_?IntegerQ, l_?IntegerQ, m_?IntegerQ, a_?NumericQ, nInv_?IntegerQ] := Module[{Alm, n},
  Alm = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]] + 2 a m \[Omega] - a^2 \[Omega]^2;
  \[Beta]freq[nInv, \[Omega], Alm, s, m, a] + ContinuedFractionK[-\[Alpha]freq[nInv-n, \[Omega], s, m, a] \[Gamma]freq[nInv-n+1, \[Omega], s, m, a], \[Beta]freq[nInv-n,\[Omega], Alm, s, m, a], {n, 1, nInv}] + CF[-\[Alpha]freq[n-1, \[Omega], s, m, a] \[Gamma]freq[n, \[Omega], s, m, a], \[Beta]freq[n, \[Omega], Alm, s, m, a], {n, nInv+1}]
];

Leaver31Ang[\[Omega]_?NumericQ, Alm_?NumericQ, s_?IntegerQ, m_?IntegerQ, a_?NumericQ, nInv_?IntegerQ] := Module[{n},
  \[Beta]ang[nInv, \[Omega], Alm, s, m, a] + ContinuedFractionK[-\[Alpha]ang[nInv-n, \[Omega], s, m, a] \[Gamma]ang[nInv-n+1, \[Omega], s, m, a], \[Beta]ang[nInv-n,\[Omega], Alm, s, m, a], {n, 1, nInv}] + CF[-\[Alpha]ang[n-1, \[Omega], s, m, a] \[Gamma]ang[n, \[Omega], s, m, a], \[Beta]ang[n, \[Omega], Alm, s, m, a], {n, nInv+1}]
];


Options[QNMFrequencyLeaver] = {"InitialGuess" -> Automatic};


QNMFrequencyLeaver[s_Integer, l_Integer, m_Integer, n_, a_, OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, inInc, nInv = n, k, funcRad, funcAng, prec = prec[a]},
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
  If[a==0,
    funcRad[\[Omega]_?NumericQ] := funcRad[\[Omega]] = Leaver[\[Omega], s, l, nInv];,
    funcRad[\[Omega]_?NumericQ] := funcRad[\[Omega]] = Leaver31[\[Omega], s, l, m, a, nInv];    
  ];
  \[Omega]QNM /. Check[FindRoot[funcRad[\[Omega]QNM]==0, {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum];
  \[Omega]QNM /. Quiet[Check[FindRoot[funcRad[\[Omega]QNM]==0, {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum], FindRoot::nlnum]
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
    "Leaver",
      \[Omega] = QNMFrequencyLeaver[s, l, m, n, a],
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyLeaver[s, l, m, n, a, opts];,
   Automatic | "SpheroidalEigenvalue",
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a];,
    {"SpheroidalEigenvalue", Rule[_,_]...},
	  opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a, opts];,
    "Spectral1D",
      If[a!=0,
        Message[QNMFrequency::nokerr, "Spectral1D", a];
        \[Omega] = $Failed;
        ,
        \[Omega] = SpectralInitialGuess[s, l, n];
      ];,
    "Spectral2D",
      If[l =!= All || n =!= All,
        Message[QNMFrequency::allln, "Spectral2D", l, n];
        \[Omega] = $Failed;,
        \[Omega] = SpectralInitialGuessKerr[s, m, a];
      ];,
    {"Spectral2D", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SpectralInitialGuessKerr]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      If[l!= All || n!= All,
        Message[QNMFrequency::allln, "Spectral2D", l, n];
        \[Omega] = $Failed;,
  	  \[Omega] = SpectralInitialGuessKerr[s, m, a, opts];
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
