load "X0p_NiceModel.m";
load "Chabauty_MWSieve.m";

//we find models for X79 and X79/w79

C := CuspForms(79);
"Dimension of CuspForms(79) is: ", Dimension(C);

AL79 := AtkinLehnerOperator(C, 79);
N79 := Nullspace(Matrix(AL79 - 1));

"Dimesion of eigenspace lambda = 1 for w79 is: ", Dimension(N79);

N79c := Nullspace(Matrix(AL79 + 1));

"Dimesion of eigenspace lambda = -1 for w79 is: ", Dimension(N79c);
"";

B79 := [&+[(Integers()!(2*Eltseq(Basis(N79)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N79)]];
B79c := [&+[(Integers()!(2*Eltseq(Basis(N79c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N79c)]];

X79 := modformeqns(B79c cat B79, 79, 500, 1);
"Nice model for X0(79) is:";
X79;
"";
RR<[u]> := CoordinateRing(AmbientSpace(X79));
n := Dimension(AmbientSpace(X79));

H := Matrix(RationalField(), 6, 6, [1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]];
w79 := iso<X79 -> X79 | rows, rows>;
"w79 on X0(79) is given by:";
w79;
"";

X79w79, quotMap := CurveQuotient(AutomorphismGroup(X79, [w79]));

"Genus of X0(79) is ", Genus(X79);
"Genus of X0(79)/w79 is ", Genus(X79w79);
"";

pts := PointSearch(X79, 2);

"We have found these points on X0(79):";
pts;
"";
 
Dtor := Divisor(pts[1]) - Divisor(pts[2]);

b, ff := IsPrincipal(13*Dtor);
"Is Dtor := pts[1] - pts[2] a generator for J0(79)(\Q)_tors? ", b and not(IsPrincipal(Dtor));
"";

E, eMap := EllipticCurve(X79w79, quotMap(pts[1]));
"X0(79)/w79 is actually the following elliptic curve:";
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
		if w79(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

"We have found ", #deg2, " points on X_0(79)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(79)/w79.";

//Finally, we do the sieve.
A := AbelianGroup([0, 13]);
divs := [D1, Dtor];
genusC := Genus(X79w79);
bp := deg2pb[1];
w79Matrix := Matrix(w79);

//proposition 3.1 - Box
I := 2;

primes:=[3,5,7,11];
B0, iA0 := sub<A | A.1>;

//leaving out representative 0 on purpose, we want to eliminate all others
W0 := {a*A.2 : a in [1..12]};

bret := MWSieve(X79, w79Matrix, genusC, primes, A, divs, I, bp, B0, iA0, W0, deg2);

"MWSieve achieved its goal? (true if succeeded, number otherwise)";
assert bret eq true;
bret;

"It follows that if there is an unknown Q in X0(79)^2(\Q), then 2[Q - bp] is fixed by w79.";
"This implies that [Q - bp] is fixed by w79 (since there is no 2-torsion)";
"Then Q ~ w79(Q), which implies that Q = w79(Q) because X0(79) is not hyperelliptic.";
"Then either Q is a pullback, or it is fixed by w79 pointwise.";
"If P = (X1:X2:X3:X4:X5:X6) is fixed by w79, then either X6 = 0 or P = (0:0:0:0:0:1) (not on X0(79))";
"";

I := IdentityMatrix(Rationals(), Genus(X79));
v := Basis(Kernel(w79Matrix + I))[1];
CR<[x]> := CoordinateRing(AmbientSpace(X79));
I := ideal<CR | &+[v[i]*x[i] : i in [1..Genus(X79)]]>;

Z := Scheme(X79, I);
"We find all such P by imposing condition X6 = 0 and finding points on the scheme:";
"";
Z;
"";

pts := PointsOverSplittingField(Z);
"All such P are:";
pts;
"";
"Clearly, none such points are quadratic.";
"Hence there are no quadratic points on X0(79) not coming from pullbacks of rationals.";