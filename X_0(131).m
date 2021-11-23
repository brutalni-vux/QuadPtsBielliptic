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

"Genus of X0(131) is ", Genus(X131);
"Genus of X0(131)/w131 is 1";
"";

P1 := X131![1,0,0,0,0,0,0,0,0,0,1];
P2 := X131![-1,0,0,0,0,0,0,0,0,0,1];

pts := [P1, P2];

"We have found these points on X0(131):";
pts;
"";
 
Dtor := Divisor(pts[1]) - Divisor(pts[2]);

b, ff := IsPrincipal(65*Dtor);
"Is Dtor := pts[1] - pts[2] a generator for J0(131)(Q)_tors? ", b and not(IsPrincipal(5*Dtor)) and not(IsPrincipal(13*Dtor));
"";

//known degree 1 places
pls1 := [Place(pts[1]), Place(pts[2])];

//known degree 2 places, we conjecture all are pullbacks
pls2:=[];

deg2:=[];
deg2pb:=[];

for i in [1..#pls1] do
	for j in [i..#pls1] do
		deg2 := Append(deg2, 1*pls1[i] + 1*pls1[j]);
		if w131(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

"We have found ", #deg2, " points on X_0(131)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(131)/w131.";

//Finally, we do the sieve.
A := AbelianGroup([65]);
divs := [Dtor];
genusC := 1;
bp := deg2pb[1];
w131Matrix := Matrix(w131);

primes := [3, 5];
B0, iA0 := sub<A | A.1>;
W0 := {0*A.1};

B, iA, W := MWSieve(X131, w131Matrix, genusC, primes, A, divs, bp, B0, iA0, W0, deg2);

"";
"For unknown Q in X0(131)^2(\Q) we have (1 - w131)[Q - bp] is in a coset of ", B, "represented by an element of ", W;
"";

"It follows that if there is an unknown Q in X0(131)^2(\Q), then (1 - w131)[Q - bp] == 0.";
"This implies that [Q - bp] is fixed by w131";
"Then Q ~ w131(Q), which implies that Q = w131(Q) because X0(131) is not hyperelliptic.";
"Then either Q is a pullback, or it is fixed by w131 pointwise.";
"If P = (X1:X2:X3:X4:X5:X6:X7:X8:X9:X10:X11) is fixed by w131,";
"then either X11 = 0 or P = (0:0:0:0:0:0:0:0:0:1) (not on X0(131))";
"";

I := IdentityMatrix(Rationals(), Genus(X131));
v := Basis(Kernel(w131Matrix + I))[1];
CR<[x]> := CoordinateRing(AmbientSpace(X131));
I := ideal<CR | &+[v[i]*x[i] : i in [1..Genus(X131)]]>;

Z := Scheme(X131, I);
"We find all such P by imposing condition X11 = 0 and finding points on the scheme:";
"";
Z;
"";

pts := PointsOverSplittingField(Z);
"All such P are:";
pts;
"";
"Clearly, none such points are quadratic.";
"Hence there are no quadratic points on X0(131) not coming from pullbacks of rationals.";