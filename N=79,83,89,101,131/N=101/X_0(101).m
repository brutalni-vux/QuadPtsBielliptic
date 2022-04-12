load "X0p_NiceModel.m";
load "Chabauty_MWSieve.m";

//we find models for X101 and X101/w101

C := CuspForms(101);
"Dimension of CuspForms(101) is: ", Dimension(C);

AL101 := AtkinLehnerOperator(C, 101);
N101 := Nullspace(Matrix(AL101 - 1));

"Dimesion of eigenspace lambda = 1 for w101 is: ", Dimension(N101);

N101c := Nullspace(Matrix(AL101 + 1));

"Dimesion of eigenspace lambda = -1 for w101 is: ", Dimension(N101c);
"";

B101 := [&+[(Integers()!(2*Eltseq(Basis(N101)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N101)]];
B101c := [&+[(Integers()!(2*Eltseq(Basis(N101c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N101c)]];

X101 := modformeqns(B101c cat B101, 101, 500, 1);
"Nice model for X0(101) is:";
X101;
"";
RR<[u]> := CoordinateRing(AmbientSpace(X101));
n := Dimension(AmbientSpace(X101));

H := Matrix(RationalField(), 8, 8, [1,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0, 0,0,1,0,0,0,0,0, 0,0,0,1,0,0,0,0, 0,0,0,0,1,0,0,0, 0,0,0,0,0,1,0,0, 0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]];
w101 := iso<X101 -> X101 | rows, rows>;
"w101 on X101 is given by:";
w101;

X101w101, quotMap := Curve_and_Map(X101, 7);

"Genus of X0(101) is ", Genus(X101);
"Genus of X0(101)/w101 is ", Genus(X101w101);
"";

P1 := X101![-1, 0, 0, 0, 0, 0, 0, 1];
P2 := X101![1, 0, 0, 0, 0, 0, 0, 1];

pts := [P1, P2];

"We have found these points on X0(101):";
pts;
"";
 
Dtor := Divisor(pts[1]) - Divisor(pts[2]);

b, ff := IsPrincipal(25*Dtor);
"Is Dtor := pts[1] - pts[2] a generator for J0(101)(Q)_tors? ", b and not(IsPrincipal(Dtor)) and not(IsPrincipal(5*Dtor));
"";

E, eMap := EllipticCurve(X101w101, quotMap(pts[1]));
"X0(101)/w101 is actually the following elliptic curve:";
E;
"";

mp := quotMap*eMap;
MWE, phi, tf1, tf2 := MordellWeilGroup(E);
"It has the following MW group:";
MWE;
"";
assert tf1 and tf2;

//we get the pullback of generator for free part of E
D := Pullback(mp, Place(phi(MWE.1)));
bp := Pullback(mp, Place(Zero(E)));
D1 := D - bp;

//known degree 1 places
pls1 := [Place(pts[1]), Place(pts[2])];

//known degree 2 places, we conjecture all are pullbacks
pls2:=[];

deg2:=[];
deg2pb:=[];

for i in [1..#pls1] do
	for j in [i..#pls1] do
		deg2 := Append(deg2, 1*pls1[i] + 1*pls1[j]);
		if w101(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

"We have found ", #deg2, " points on X_0(101)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(101)/w101.";

//Finally, we do the sieve.
A := AbelianGroup([0, 25]);
divs := [D1, Dtor];
genusC := Genus(X101w101);
bp := deg2pb[1];
w101Matrix := Matrix(w101);

//proposition 3.1 - Box
I := 2;

primes:=[3,5,7,11,13,17,19,23];
B0, iA0 := sub<A | A.1>;

//leaving out representative 0 on purpose, we want to eliminate all others
W0 := {a*A.2 : a in [1..24]};

bret := MWSieve(X101, w101Matrix, genusC, primes, A, divs, I, bp, B0, iA0, W0, deg2);

"MWSieve achieved its goal? (true if succeeded, number otherwise)";
assert bret eq true;
bret;

"It follows that if there is an unknown Q in X0(101)^2(\Q), then 2[Q - bp] is fixed by w101.";
"This implies that [Q - bp] is fixed by w101 (since there is no 2-torsion)";
"Then Q ~ w101(Q), which implies that Q = w101(Q) because X0(101) is not hyperelliptic.";
"Then either Q is a pullback, or it is fixed by w101 pointwise.";
"If P = (X1:X2:X3:X4:X5:X6:X7:X8) is fixed by w101, then either X8 = 0 or P = (0:0:0:0:0:0:0:1) (not on X0(101))";
"";

I := IdentityMatrix(Rationals(), Genus(X101));
v := Basis(Kernel(w101Matrix + I))[1];
CR<[x]> := CoordinateRing(AmbientSpace(X101));
I := ideal<CR | &+[v[i]*x[i] : i in [1..Genus(X101)]]>;

Z := Scheme(X101, I);
"We find all such P by imposing condition X8 = 0 and finding points on the scheme:";
"";
Z;
"";

pts := PointsOverSplittingField(Z);
"All such P are:";
pts;
"";
"Clearly, none such points are quadratic.";
"Hence there are no quadratic points on X0(101) not coming from pullbacks of rationals.";
