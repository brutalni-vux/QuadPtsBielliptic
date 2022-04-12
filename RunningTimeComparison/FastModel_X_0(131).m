load "X0p_NiceModel.m";
load "Chabauty_MWSieve_131.m";

//we find models for X131 and X131/w131

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

"";
RR<[u]> := CoordinateRing(AmbientSpace(X131));
n := Dimension(AmbientSpace(X131));

H := Matrix(RationalField(), 11, 11, [1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]] ;
w131 := iso<X131 -> X131 | rows, rows>;
"w131 on X131 is given by:";
w131;
"";


X131w131, quotMap := Curve_and_Map(X131, 10);

P1 := X131![1,0,0,0,0,0,0,0,0,0,1];
P2 := X131![-1,0,0,0,0,0,0,0,0,0,1];

pts := [P1, P2];

E, eMap := EllipticCurve(X131w131, quotMap(pts[1]));
"X0(131)/w131 is actually the following elliptic curve:";
E;
"";

mp := quotMap*eMap;

"Quotient map is:";
mp;
