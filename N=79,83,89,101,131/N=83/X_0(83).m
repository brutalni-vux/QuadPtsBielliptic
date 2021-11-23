load "X0p_NiceModel.m";
load "Chabauty_MWSieve.m";

//we find models for X83 and X83/w83

C := CuspForms(83);
"Dimension of CuspForms(83) is: ", Dimension(C);

AL83 := AtkinLehnerOperator(C, 83);
N83 := Nullspace(Matrix(AL83 - 1));

"Dimesion of eigenspace lambda = 1 for w83 is: ", Dimension(N83);

N83c := Nullspace(Matrix(AL83 + 1));

"Dimesion of eigenspace lambda = -1 for w83 is: ", Dimension(N83c);
"";

B83 := [&+[(Integers()!(2*Eltseq(Basis(N83)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N83)]];
B83c := [&+[(Integers()!(2*Eltseq(Basis(N83c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N83c)]];

X83 := modformeqns(B83c cat B83, 83, 500, 1);
"Nice model for X0(83) is:";
X83;
"";
RR<[u]> := CoordinateRing(AmbientSpace(X83));
n := Dimension(AmbientSpace(X83));

H := Matrix(RationalField(), 7, 7, [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]] ;
w83 := iso<X83 -> X83 | rows, rows>;
"w83 on X0(83) is given by:";
w83;
"";

X83w83, quotMap := CurveQuotient(AutomorphismGroup(X83, [w83]));

"Genus of X0(83) is ", Genus(X83);
"Genus of X0(83)/w83 is ", Genus(X83w83);
"";

pts := PointSearch(X83, 2);

"We have found these points on X0(83):";
pts;
"";
 
Dtor := Divisor(pts[1]) - Divisor(pts[2]);

b, ff := IsPrincipal(41*Dtor);
"Is Dtor := pts[1] - pts[2] a generator for J0(83)(Q)_tors? ", b and not(IsPrincipal(Dtor));
"";

E, eMap := EllipticCurve(X83w83, quotMap(pts[1]));
"X0(83)/w83 is actually the following elliptic curve:";
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
		if w83(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

"We have found ", #deg2, " points on X_0(83)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(83)/w83.";

//Finally, we do the sieve.
A := AbelianGroup([0, 41]);
divs := [D1, Dtor];
genusC := Genus(X83w83);
bp := deg2pb[1];
w83Matrix := Matrix(w83);

//proposition 3.1 - Box
I := 2;

primes:=[3,5,7,11,13,17,19,23];
B0, iA0 := sub<A | A.1>;

//leaving out representative 0 on purpose, we want to eliminate all others
W0 := {a*A.2 : a in [1..40]};

bret := MWSieve(X83, w83Matrix, genusC, primes, A, divs, I, bp, B0, iA0, W0, deg2);

"MWSieve achieved its goal? (true if succeeded, number otherwise)";
assert bret eq true;
bret;

"It follows that if there is an unknown Q in X0(83)^2(\Q), then 2[Q - bp] is fixed by w83.";
"This implies that [Q - bp] is fixed by w83 (since there is no 2-torsion)";
"Then Q ~ w83(Q), which implies that Q = w83(Q) because X0(83) is not hyperelliptic.";
"Then either Q is a pullback, or it is fixed by w83 pointwise.";
"If P = (X1:X2:X3:X4:X5:X6:X7) is fixed by w83, then either X7 = 0 or P = (0:0:0:0:0:0:1) (not on X0(83))";
"";

I := IdentityMatrix(Rationals(), Genus(X83));
v := Basis(Kernel(w83Matrix + I))[1];
CR<[x]> := CoordinateRing(AmbientSpace(X83));
I := ideal<CR | &+[v[i]*x[i] : i in [1..Genus(X83)]]>;

Z := Scheme(X83, I);
"We find all such P by imposing condition X7 = 0 and finding points on the scheme:";
"";
Z;
"";

pts := PointsOverSplittingField(Z);
"All such P are:";
pts;
"";
"Clearly, none such points are quadratic.";
"Hence there are no quadratic points on X0(83) not coming from pullbacks of rationals.";
