load "X0p_NiceModel.m";
load "Chabauty_MWSieve_131.m";

//compute model for X_0(131) using the diagonal cuspform basis (such that the Atkin-Lehner matrix is diagonal)

C := CuspForms(131);
"Dimension of CuspForms(131) is: ", Dimension(C);

AL131 := AtkinLehnerOperator(C, 131);
N131 := Nullspace(Matrix(AL131 - 1));

"Dimesion of eigenspace lambda = 1 for w131 is: ", Dimension(N131);

N131c := Nullspace(Matrix(AL131 + 1));

"Dimesion of eigenspace lambda = -1 for w131 is: ", Dimension(N131c);
"";

B131 := [&+[(Integers()!(2*Eltseq(Basis(N131)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N131)]];
B131c := [&+[(Integers()!(2*Eltseq(Basis(N131c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N131c)]];

X131 := modformeqns(B131c cat B131, 131, 500, 1);
"Nice model for X0(131) is:";
X131;
