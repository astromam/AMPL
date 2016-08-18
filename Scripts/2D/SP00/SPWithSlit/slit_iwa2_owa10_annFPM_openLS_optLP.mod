# square aperture, bounded box focal plane.
# with apodizer array as free parameter,
# maximize apodizer transmission subject to Lyot plane field constraint

param pi := 4*atan(1);

param IWA := 2.0;
param OWA := 10;
#param s := 3.;
param s := 4.2423; # linear prolate eigenvalue for FPM spot radius 2 lam/D

param M_p := 1000;			# discretization (first pupil plane)
param M_q := 1000;			# discretization (Lyot plane)
param M_xi := 8*OWA;		# discretization (first focal plane)

param p {j in 1..M_p} := 0.5*(j-0.5)/M_p;
param dp := p[2]-p[1];
param q {j in 1..M_q} := 0.5*(j-0.5)/M_q;
param dq := q[2]-q[1];

var A_SP {j in 1..M_p} <=1, >=0, binary, :=0.5;

set Xis ordered;
let Xis := setof {j in 1..M_xi} OWA*(j-0.5)/M_xi;
param dxi := OWA/M_xi;
set fp1 := setof {xi in Xis: xi>=IWA} xi;

var E_fp1 {xi in fp1} = sum {j in 1..M_p} 2*A_SP[j]*cos(2*pi*xi*p[j])*dp;
var E_Lp {j in 1..M_q} = sum {xi in fp1} 2*E_fp1[xi]*cos(2*pi*q[j]*xi)*dxi;
var E_fp1_peak = 2 * sum {j in 1..M_p} A_SP[j]*dp;

maximize thrupt: 2 * sum {j in 1..M_p} A_SP[j]*dp;
subject to Lyot_constr_pos {j in 1..M_q}: E_Lp[j] <= 10^(-s)*E_fp1_peak;
subject to Lyot_constr_neg {j in 1..M_q}: -10^(-s)*E_fp1_peak <= E_Lp[j];

option solver loqo;
option loqo_options "verbose=2";
#option solver gurobix;
#option gurobi_options "outlev=1 presolve=0 method=1 crossover=0";

solve;

printf {j in 1..M_p}: "%15e %15e \n", p[j], A_SP[j] > "SP_slit_iwa2_owa10_annFPM_openLS_optLP42.dat";
