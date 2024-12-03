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


ClearAttributes[{QNMFrequency, QNMFunction}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, n, a] computes the quasinormal mode frequncy";
QNMFunction::usage = "QNMFunction[\[Omega],s,m,a] computes the radial eigenfunction associated with the eigenvalue \[Omega]"; 


(* ::Subsection:: *)
(*Error messages*)


QNMFrequency::real = "Only real values of a are allowed, but a=`1`.";


(* ::Subsection:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*Calculate QNM Frequency*)


(*CP comment: To do: currently just n=0. Generalise?*)


(* ::Subsection::Closed:: *)
(*Parameters*)


(*CP: cleared this up*)
M=1;(*re-instate at some point????*)


(* ::Subsection::Closed:: *)
(*Operator Functions*)


\[CapitalDelta][\[Rho]_,a_]:=1-2M \[Rho]+a^2 \[Rho]^2;
A[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=2 I \[Omega]-2(1+s)\[Rho]+2(I \[Omega](a^2-8M^2)+I m a +(s+3)M)\[Rho]^2+4(2 I \[Omega] M-1) a^2 \[Rho]^3;
B[\[Omega]_,s_Integer,m_Integer,a_,\[Rho]_]:=(a^2-16M^2)\[Omega]^2+2(m a +2 I s M)\[Omega]+2(4(a^2-4M^2)M \[Omega]^2+(4 m a M-4I (s+2) M^2+I a^2)\[Omega]+ I m a+(s+1)M)\[Rho]+2(8 M^2 \[Omega]^2+6 I M \[Omega]-1)a^2 \[Rho]^2;
M\[Rho][\[Omega]_,s_Integer,m_Integer,a_]:=-\[Rho]^2\[CapitalDelta][\[Rho],a]R''[\[Rho]]+A[\[Omega],s,m ,a,\[Rho]]R'[\[Rho]]+B[\[Omega],s,m ,a,\[Rho]]R[\[Rho]];


(* ::Subsection::Closed:: *)
(*Discretization*)


rmax[a_]:=If[a==0, 1/2,1/a^2 (1-Sqrt[1-a^2]) ];
(* Define the grid *)
\[Rho]grid[a_]:=rmax[a] Sort[N[1/2 (1+Cos[\[Pi] Range[0,1,1/(Nz-1)]])]];
(* Get differentiation matrices based on the grid *)
Dgrid[\[Rho]grid_]:= Dgrid[\[Rho]grid]= NDSolve`FiniteDifferenceDerivative[Derivative[1],{\[Rho]grid},
	DifferenceOrder->{"Pseudospectral"},PeriodicInterpolation->{False}]["DifferentiationMatrix"];
DDgrid[\[Rho]grid_]:= DDgrid[\[Rho]grid] = NDSolve`FiniteDifferenceDerivative[Derivative[2],{\[Rho]grid},
	DifferenceOrder->{"Pseudospectral"},PeriodicInterpolation->{False}]["DifferentiationMatrix"];


(* ::Subsection::Closed:: *)
(*Radial Eigenvalues*)


(* Compute the eigenvalue of the Radial equation *)
\[Lambda]r[\[Omega]guess_,s_Integer, m_Integer, a_]:=Module[
	{\[Rho]grida,dd1,dd2,Mat,DiscretizationRules},
		(* Call grids *)
		\[Rho]grida = \[Rho]grid[a];
		dd1 = Dgrid[\[Rho]grida];
		dd2 = DDgrid[\[Rho]grida];
		(* Define discretization in terms of grids *)
		DiscretizationRules={\[Rho]->\[Rho]grida,R''[\[Rho]]->dd2, 
			R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[Nz]};
		(* Evaluate discretized operator *)
		Mat=M\[Rho][\[Omega]guess,s,m,a]/.DiscretizationRules;
		(* Return sorted eigenvalues *)
		Eigenvalues[{Mat}]//SortBy[#,Norm]&
]


(* CP comment: Computes the difference between the  radial and angular separation constant! *)
\[Delta]\[Lambda][\[Omega]guess_,s_Integer,l_Integer,m_Integer,a_]:=Module[
	{\[CapitalLambda]},
		\[CapitalLambda]=SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
		-\[Lambda]r[\[Omega]guess,s,m,a][[1]]+(\[CapitalLambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2)
]


(* ::Subsection::Closed:: *)
(*Newton Raphson solver*)


(*CP commment: what the algorithm solves for is not the SpinWeightedSpheroidalEigenvalue, but the combination (SpinWeightedSpheroidalEigenvalue+2 m a \[Omega]guess- (a \[Omega]guess)^2) *)
(*CP: what do we want to output here?*)
QNMFrequency[s_Integer,l_Integer,m_Integer,a_,tol_:10^-6]:=Module[
	{\[Lambda],\[Omega]guess,F,Fp,\[Epsilon],\[Gamma],\[Delta]\[Omega]},
		(* Check for real spin *)
		If[!RealValuedNumberQ[a],
		Message[QNMFrequency::real, a];
		Return[$Failed];
		];
  
		(* Intial guess *)
		\[Omega]guess = 1;
		\[Delta]\[Omega]=1;
		\[Epsilon]=10^-8;(*better way?*)
		\[Gamma]=1;(*allow this to change?*)
		
		(* Iterate until tolerance is reached *)
		Until[Norm[\[Delta]\[Omega]]<tol,
			F=\[Delta]\[Lambda][\[Omega]guess, s,l, m, a];
			Fp=(\[Delta]\[Lambda][\[Omega]guess+\[Epsilon], s,l, m, a]-\[Delta]\[Lambda][\[Omega]guess-\[Epsilon], s,l, m, a])/(2 \[Epsilon])+(\[Delta]\[Lambda][\[Omega]guess+I \[Epsilon], s,l, m, a]-\[Delta]\[Lambda][\[Omega]guess-I \[Epsilon], s,l, m, a])/(2 \[Epsilon] I);
			\[Omega]guess= \[Omega]guess-(\[Gamma] F)/Fp;
			\[Delta]\[Omega]= -(\[Gamma] F)/Fp
		];
		
		(* Once result is obtained, return eigenvalue & separation constant *)
		\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
		<|"QNMfrequency"->\[Omega]guess,"SeparationConstant"->\[Lambda]|>
]


(* ::Section:: *)
(*Calculate radial function*)


(* ::Subsection:: *)
(*Use the eigenvalue to calculate the eigenfunction*)


(*CP comment: this is error prone: it takes in the QNM frequency, which the uses might give erroneously.. Maybe the two function are put together into 1, or the second one expands in content?*)
(*CP comment: There should be an option: `Hyperboloidal' or `BL'*)
QNMFunction[\[Omega]_,s_Integer,m_Integer,a_]:=Module[
	{ef,\[Rho]grida,dd1,dd2,DiscretizationRules,Mat},
		(* Call the discretized system *)
		\[Rho]grida = \[Rho]grid[a];
		dd1 = Dgrid[\[Rho]grida];
		dd2 = DDgrid[\[Rho]grida];
		(* Define discretization in terms of grids *)
		DiscretizationRules={\[Rho]->\[Rho]grida,R''[\[Rho]]->dd2, 
			R'[\[Rho]]->dd1, R[\[Rho]]->IdentityMatrix[Nz]};
		(* Evaluate discretized operator *)
		Mat=M\[Rho][\[Omega],s,m,a]/.DiscretizationRules;
		(* Calculate the eigensystem given the eigenvalue *)
		ef=Sort[Transpose[Eigensystem[{Mat}]], Norm[#1[[1]]]<Norm[#2[[1]]]&][[1,2]];
		(* Return table of radial gridpoints and eigenfunction values *)
		Table[{\[Rho]grida[[i]], ef[[i]]},{i,1,Length[\[Rho]grida]}]
]


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency, QNMFunction}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
