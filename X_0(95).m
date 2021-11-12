load "quadptssieve.m";
load "X0N_NiceModel.m";

//we find models for X0(95) and X0(95)/w19

C:=CuspForms(95);
"Dimension of CuspForms(95) is:";
Dimension(C);

w19:=AtkinLehnerOperator(C, 19);
N19:=Nullspace(Matrix(w19 - 1));

"Dimesion of eigenspace lambda=1 for w19 is:";
Dimension(N19);

N19c:=Nullspace(Matrix(w19 + 1));

"Dimesion of eigenspace lambda=-1 for w19 is:";
Dimension(N19c);

B19:=[&+[(Integers()!(2*Eltseq(Basis(N19)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N19)]];
B19c:=[&+[(Integers()!(2*Eltseq(Basis(N19c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N19c)]];

X95,jMap:= modformeqns(B19c cat B19, 95, 500, 19);
"Model for X95 is:";
X95;
"";
RR<[u]>:=CoordinateRing(AmbientSpace(X95));
n:=Dimension(AmbientSpace(X95));

H:=Matrix(RationalField(), 9, 9, [1,0,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0, 0,0,1,0,0,0,0,0,0, 0,0,0,1,0,0,0,0,0, 0,0,0,0,1,0,0,0,0, 0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,-1,0,0, 0,0,0,0,0,0,0,-1,0, 0,0,0,0,0,0,0,0,-1 ]);
rows:=[[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]] ;
w19:=iso<X95->X95 | rows,rows>;
"w19 on X0(95) is given by:";
w19;
"";

X95w19, quotMap:=CurveQuotient(AutomorphismGroup(X95,[w19]));
"Genus of X0(95) is ", Genus(X95);
"Genus of X0(95)/w19 is ", Genus(X95w19);
"";

//some pts on X095
K19<r19>:=QuadraticField(-19);

P1:=X95![-1,0,0,0,0,0,1,0,0];
P2:=X95![1,0,0,0,0,0,1,0,0];
P3:=X95![-3/5,0,-2/5,-1/5,-1/5,-1/5,1,0,0];
P4:=X95![3/5,0,2/5,1/5,1/5,1/5,1,0,0];

"By finding poles of j-map, we find that we have these 4 cusps:";
Poles(jMap);
"";

P5:=X95(K19)![1/14*(-r19-17), 1/7*(r19-11), 1/7*(r19+3), 1/7*(r19-4), 1/14*(r19+3),1,0,0,0];
P6:=X95(K19)![1/14*(r19-17), 1/7*(-r19-11), 1/7*(-r19+3), 1/7*(-r19-4), 1/14*(-r19+3),1,0,0,0];

"One quadratic but not rational point is P5:", P5;
"Is P5 fixed by w19?";
w19(P5) eq P5;
"Hence, P5 can't come from a pullback of rational point.";
"";

D1:=Divisor(P1) - Divisor(P3);
bb, ff := IsPrincipal(180*D1);
"Is D1:=Divisor(P1) - Divisor(P3) of order 180?";
bb and not (IsPrincipal(90*D1) or IsPrincipal(60*D1) or IsPrincipal(36*D1));
"";

D2:=Divisor(P1)-Divisor(P2)-10*D1;
"Is D2:=Divisor(P1)-Divisor(P2)-10*D1 of order 6?";
bb, ff := IsPrincipal(6*D2);
bb and not(IsPrincipal(3*D2) or IsPrincipal(2*D2));
"";

//since D2 is of order 6, its D1-part is 30*k*D1 for k in 0 to 5
//3*D2 has D1-part equal to 0 or 90*D1
//2*D2 has D1-part equal to 0, 60*D1 or -60*D1
//we check that D1-free part of D2 is of order 6
"Is D2:=Divisor(P1)-Divisor(P2)-10*D1 of order 6 in J(Q)_tor/<D1>?";
not(IsPrincipal(3*D2) or IsPrincipal(3*D2 - 90*D1) or IsPrincipal(2*D2) or IsPrincipal(2*D2-60*D1) or IsPrincipal(2*D2 + 60*D1));
"Hence, D2 and D1 generate a subgroup isomorphic to Z/6Z + Z/180Z";
"";

X3:=ChangeRing(X95,GF(3));
CG,phi,psi:=ClassGroup(X3);
Z:=FreeAbelianGroup(1);
degr:=hom<CG->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(CG)]>;
JF3:=Kernel(degr); // This is isomorphic to J_X(\F_3).
"J(Q)_tor is isomorphic to a subgroup of:";
JF3;
"";

X7:=ChangeRing(X95,GF(7));
CG,phi,psi:=ClassGroup(X7);
Z:=FreeAbelianGroup(1);
degr:=hom<CG->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(CG)]>;
JF7:=Kernel(degr); // This is isomorphic to J_X(\F_7).
"J(Q)_tor is isomorphic to a subgroup of:";
JF7;
"Hence, J(Q)_tor doesn't have two independent elements of order 4.";
"";

X11:=ChangeRing(X95,GF(11));
CG,phi,psi:=ClassGroup(X11);
Z:=FreeAbelianGroup(1);
degr:=hom<CG->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(CG)]>;
JF11:=Kernel(degr); // This is isomorphic to J_X(\F_11).
"J(Q)_tor is isomorphic to a subgroup of:";
JF11;
"Hence, J(Q)_tor doesn't have two independent elements of order 10.";
"";

"Since we have found a torsion subgroup iso. to Z/6Z + Z/180Z, and rk(J0(95)(Q)) == 0,";
"we know that J(Q) is a subgroup of:";
"Z/2Z + Z/2Z + Z/6Z + Z/180Z";
"";

//known degree 1 places
pls1:=[Place(P1), Place(P2), Place(P3), Place(P4)];

//known degree 2 places, we know P5 is not a pullback of a rational
pls2:=[Place(P5)];

deg2 := [1*pl : pl in pls2];
deg2pb:=[];

for i in [1..#pls1] do
	for j in [i..#pls1] do
		deg2:=Append(deg2, 1*pls1[i] + 1*pls1[j]);
		if w19(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb:=Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;


/*deg2:=[1*pl : pl in pls2] cat [1*pl1 + 1*pl2 : pl1 in pls1, pl2 in pls1];
deg2pb:=[1*pl1 + 1*pl2 : pl1 in pls1, pl2 in pls1 | w19(RepresentativePoint(pl1)) eq RepresentativePoint(pl2)];*/

deg2npb:=[DD : DD in deg2 | not DD in deg2pb];

"We have found ", #deg2, " points on X_0(95)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(95)/w19.";
#deg2npb, "of them are non-pullbacks";

//Finally, we do the sieve.
A:=AbelianGroup([6, 180]);
divs:=[D2, D1];
genusC:=Genus(X95w19);
auts:=[H];
bp:=deg2pb[1];

//We have shown that 2J(Q)_tor <= A
I:=2;

primes := [11, 13];

"Succeeded in proving that we have all exceptional quadratic pts? (true if succeeded, number otherwise)";
bret := MWSieve(deg2,primes,X95,A,divs,auts,genusC,deg2pb,deg2npb,I,bp);
bret;
assert bret eq true;
"";

"Hence, there are no quadratic points on X0(95) not coming from X0(95)/w19(Q).";
"We now have to find rational points on X0(95)/w19(Q) and check their pullbacks.";
"";

Hyp, mp := SimplifiedModel(X95w19);
"X0(95)/w19 is actually ", Hyp;

"Rank of J(X0(95)/w19)(Q) is", RankBound(Jacobian(Hyp));
"Hence rk(J(X0(95)/w19)(Q)) = rk(J0(95)(Q)) = 0.";
"";

ptsMax := #Jacobian(ChangeRing(Hyp, GF(3)));
"J(X0(95)/w19)(Q)_tors has at most ", ptsMax, " points.";

K5<w> := QuadraticField(5);
Hyp5:=ChangeRing(Hyp, K5);

Inf1 := PointsAtInfinity(Hyp5)[1];
Inf2 := PointsAtInfinity(Hyp5)[2];
Quad1 := Hyp5![1/2*(-w + 3), 1/2*(5*w - 7)];
Quad2 := Hyp5![1/2*(w + 3), 1/2*(-5*w - 7)];

Div1 := Inf1 - Inf2;
Div2 := (Quad1 - Inf1) + (Quad2 - Inf2);

S := {};

for i := 1 to 10 do
	for j := 1 to 10 do
		S := S join {i*Div1 + j*Div2};
	end for;
end for;

"We have found ", #S, " points on J(X0(95)/w19)(Q) out of possible ", ptsMax;
assert #S eq ptsMax;
"Hence, we know whole J(X0(95)/w19)(Q).";
"";

R:={};
for a in S do
	R:=R join SequenceToSet(Roots(a[1]));
end for;

"Using Mumford representations, we get that x coordinates of all rational";
"plus non-obvious quadratic points on X0(95)/w19 are (excluding pts at infinity):";

for rt in R do
	rt[1];
end for;

"Hence, only rational points are points at infinity";
"";
for i in [1..#PointsAtInfinity(Hyp)] do
	"Pullback of point ", PointsAtInfinity(Hyp)[i], " is:";
	Decomposition(Pullback(quotMap*mp, Place(PointsAtInfinity(Hyp)[i])));
end for;

"Hence, the only quadratic point up to Galois conjugacy, apart from rational cusps, is ", P5;