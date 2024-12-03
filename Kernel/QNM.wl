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


ClearAttributes[{QNMFrequency, QNMMode}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, a, tol][NN] computes the quasinormal mode frequncy following arXiv:2202.03837 up to a tolerance tol with NN gridpoints";
QNMMode::usage = "QNMMode[s,l,m,a,coords,tol][NN] computes the frequency \[Omega] and radial eigenfunction, in either Hyperboloidal or BL coordinates"; 


(* ::Subsection:: *)
(*Error messages*)


QNMFrequency::real = "Only real values of a are allowed, but a=`1`.";
QNMMode::coords = "Must specify coordinates as either 'Hyperboloidal' or 'BL'";


(* ::Subsection:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*Calculate QNM Frequency*)


(*CP comment: To do: currently just n=0. Generalise?*)


(* ::Subsection:: *)
(*Parameters*)


(*CP: cleaned this up*)
M=1;(*re-instate at some point????*)


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
		-\[Lambda]r[\[Omega]guess,s,m,a,NN][[1]]+(\[CapitalLambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2)
]


(* ::Subsection:: *)
(*Newton Raphson solver*)


(* Overloads *)
QNMFrequency[s_Integer,l_Integer,m_Integer,a_,tol_:10^-6][NN_Integer:100]:= QNMFrequency[s,l,m,a,tol,100]
QNMFrequency[s_Integer,l_Integer,m_Integer,a_,tol_:10^-6]:= QNMFrequency[s,l,m,a,tol,100]
QNMFrequency[s_Integer,l_Integer,m_Integer,a_][NN_Integer:100]:=QNMFrequency[s,l,m,a,10^-6,NN]
QNMFrequency[s_Integer,l_Integer,m_Integer,a_]:=QNMFrequency[s,l,m,a,10^-6,100]


(*CP commment: the algorithm solves for the separation constant, while SpinWeightedSpheroidalEigenvalue gives the Angular eigenvalue: SeparationConstant= (SpinWeightedSpheroidalEigenvalue+2 m a \[Omega]guess- (a \[Omega]guess)^2) *)
QNMFrequency[s_Integer,l_Integer,m_Integer,a_,tol_:10^-6,NN_Integer:100]:=Module[
	{\[Lambda],\[Omega]guess,F,Fp,\[Epsilon],\[Gamma],\[Delta]\[Omega]},
		(* Check for real spin *)
		If[!RealValuedNumberQ[a],
		Message[QNMFrequency::real, a];
		Return[$Failed];
		];
		(* Check for positive number of gridpoints *)
		If[NN <= 0, NN=100, Nothing];
		(* Intial guess *)
		
		\[Omega]guess = 1;
		\[Delta]\[Omega]=1;
		\[Epsilon]=10^-8;(*better way?*)
		\[Gamma]=1;(*allow this to change?*)
		
		(* Iterate until tolerance is reached *)
		Until[Norm[\[Delta]\[Omega]]<tol,
			F=\[Delta]\[Lambda][\[Omega]guess, s,l, m, a,NN];
			Fp=(\[Delta]\[Lambda][\[Omega]guess+\[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-\[Epsilon], s,l, m, a,NN])/(2 \[Epsilon])+(\[Delta]\[Lambda][\[Omega]guess+I \[Epsilon], s,l, m, a,NN]-\[Delta]\[Lambda][\[Omega]guess-I \[Epsilon], s,l, m, a,NN])/(2 \[Epsilon] I);
			\[Omega]guess= \[Omega]guess-(\[Gamma] F)/Fp;
			\[Delta]\[Omega]= -(\[Gamma] F)/Fp
		];
		
		(* Once result is obtained, return eigenvalue & separation constant *)
		\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]guess];
		<|"QNMfrequency"->\[Omega]guess,"SeparationConstant"->\[Lambda]+2 m a \[Omega]guess- (a \[Omega]guess)^2|>
]


(* ::Section:: *)
(*Calculate radial function*)


(* ::Subsection:: *)
(*Use the eigenvalue to calculate the eigenfunction*)


(*CP comment: this is error prone: it takes in the QNM frequency, which the user might give erroneously..*)
(*CP comment: conventions fixed \[CapitalDelta]^-s/r, \[Rho]=1/r see eq 6 in 2202.03837*)
QNMFunction[\[Omega]_,s_Integer,m_Integer,a_][NN_Integer]:=Module[
	{ef,\[Rho]grida,dd1,dd2,DiscretizationRules,Mat, Delta},
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
		ef=Sort[Transpose[Eigensystem[{Mat}]], Norm[#1[[1]]]<Norm[#2[[1]]]&][[1,2]];
		Delta=Table[(\[Rho]grida[[i]]^2 a^2-2M \[Rho]grida[[i]]+1)/\[Rho]grida[[i]]^2 ,{i,1,Length[\[Rho]grida]}];
		(* Return table of radial gridpoints and eigenfunction values *)
		Table[{\[Rho]grida[[i]], Delta[[i]]^(-s) \[Rho]grida[[i]]  ef[[i]]},{i,1,Length[\[Rho]grida]}]
]

(*CP comment: this is a way around the problem above. For `BL' multiply with E^(-I \[Omega] h), h=-r-4M Log r= -1/\[Rho]+4M Log \[Rho]*)
QNMMode[s_Integer,l_integer,m_Integer,a_,coords_,tol_:10^-6][NN_Integer:100]:=Module[
	{\[Omega], ef,\[Rho]grida,dd1,dd2,DiscretizationRules,Mat,RadialProfile, Delta,h},
		If[!StringQ[coords],
		Message[QNMMode::coords];
		Return[$Failed];
		];
		(* Get the frequency *)
		\[Omega]=QNMFrequency[s,l,m,a,tol][NN];
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
		ef=Sort[Transpose[Eigensystem[{Mat}]], Norm[#1[[1]]]<Norm[#2[[1]]]&][[1,2]];
		(* Return table of radial gridpoints and eigenfunction values *)
		Delta=Table[(\[Rho]grida[[i]]^2 a^2-2M \[Rho]grida[[i]]+1)/\[Rho]grida[[i]]^2 ,{i,1,Length[\[Rho]grida]}];
		h=Table[-1/\[Rho]grida[[i]]+4 M Log[\[Rho]grida[[i]]]  ,{i,1,Length[\[Rho]grida]}];
		
		If[coords=="Hyperboloidal", RadialProfile=Table[{\[Rho]grida[[i]],Delta[[i]]^-s \[Rho]grida[[i]] ef[[i]]},{i,1,Length[\[Rho]grida]}]];
		If[coords=="BL", RadialProfile=Table[{\[Rho]grida[[i]],Delta[[i]]^-s \[Rho]grida[[i]] Exp[-I \[Omega] h[[i]]] ef[[i]]},{i,1,Length[\[Rho]grida]}]];
		<|"QNMfrequency"->\[Omega],"QNMRadialProfile"->RadialProfile|>
]


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency, QNMMode}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
