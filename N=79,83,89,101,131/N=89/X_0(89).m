load "X0p_NiceModel.m";
load "Chabauty_MWSieve.m";

//we find models for X89 and X89/w89

C := CuspForms(89);
"Dimension of CuspForms(89) is: ", Dimension(C);

AL89 := AtkinLehnerOperator(C, 89);
N89 := Nullspace(Matrix(AL89 - 1));

"Dimesion of eigenspace lambda = 1 for w89 is: ", Dimension(N89);

N89c := Nullspace(Matrix(AL89 + 1));

"Dimesion of eigenspace lambda = -1 for w89 is: ", Dimension(N89c);
"";

B89 := [&+[(Integers()!(2*Eltseq(Basis(N89)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N89)]];
B89c := [&+[(Integers()!(2*Eltseq(Basis(N89c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N89c)]];

X89 := modformeqns(B89c cat B89, 89, 500, 1);
"Nice model for X0(89) is:";
X89;
"";
RR<[u]> := CoordinateRing(AmbientSpace(X89));
n := Dimension(AmbientSpace(X89));

H := Matrix(RationalField(), 7, 7, [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]] ;
w89 := iso<X89 -> X89 | rows, rows>;
"w89 on X0(89) is given by:";
w89;
"";

X89w89, quotMap := CurveQuotient(AutomorphismGroup(X89, [w89]));

"Genus of X0(89) is ", Genus(X89);
"Genus of X0(89)/w89 is ", Genus(X89w89);
"";

pts := PointSearch(X89, 2);

"We have found these points on X0(89):";
pts;
"";
 
Dtor := Divisor(pts[1]) - Divisor(pts[2]);

b, ff := IsPrincipal(22*Dtor);
"Is Dtor := pts[1] - pts[2] a generator for J0(89)(Q)_tors? ", b and not(IsPrincipal(2*Dtor)) and not(IsPrincipal(11*Dtor));
"";

E, eMap := EllipticCurve(X89w89, quotMap(pts[1]));
"X0(89)/w89 is actually the following elliptic curve:";
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
		if w89(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

"We have found ", #deg2, " points on X_0(89)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(89)/w89.";

//Finally, we do the sieve.
A := AbelianGroup([0, 22]);
divs := [D1, Dtor];
genusC := Genus(X89w89);
bp := deg2pb[1];
w89Matrix := Matrix(w89);

//proposition 3.1 - Box
I := 2;

primes := [3,5,7];
B0, iA0 := sub<A | A.1>;
W0 := {a*A.2 : a in [0..21]};

B, iA, W := MWSieve(X89, w89Matrix, genusC, primes, A, divs, I, bp, B0, iA0, W0, deg2);

//the goal is to get that 2[Q - bp] == 2k*D1, if sieve returned B and W such that
//B and W both consist of elements of the form 2kD1, that goal will be proved
"";
"Is it true that for unknown Q in X0(89)^2(\Q) we have 2[Q - bp] == 2k*D1?";
(B subset sub<A | 2*A.1>) and (Seqset(W) subset sub<A | 2*A.1>);
"";

"It follows that if there is an unknown Q in X0(89)^2(\Q), then 2[Q - bp] is fixed by w89.";
"This implies that [Q - bp] is fixed by w89 (since 2[Q - bp] = 2k*D1, see article)";
"Then Q ~ w89(Q), which implies that Q = w89(Q) because X0(89) is not hyperelliptic.";
"Then either Q is a pullback, or it is fixed by w89 pointwise.";
"If P = (X1:X2:X3:X4:X5:X6:X7) is fixed by w89, then either X7 = 0 or P = (0:0:0:0:0:0:1) (not on X0(89))";
"";

I := IdentityMatrix(Rationals(), Genus(X89));
v := Basis(Kernel(w89Matrix + I))[1];
CR<[x]> := CoordinateRing(AmbientSpace(X89));
I := ideal<CR | &+[v[i]*x[i] : i in [1..Genus(X89)]]>;

Z := Scheme(X89, I);
"We find all such P by imposing condition X7 = 0 and finding points on the scheme:";
"";
Z;
"";

pts := PointsOverSplittingField(Z);
"All such P are:";
pts;
"";
"Clearly, none such points are quadratic.";
"Hence there are no quadratic points on X0(89) not coming from pullbacks of rationals.";
