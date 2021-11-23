load "quadptssieve.m";
load "X0N_NiceModel.m";

//we find models for X119 and X119/w17

C:=CuspForms(119);
"Dimension of CuspForms(119) is:";
Dimension(C);

AL17:=AtkinLehnerOperator(C, 17);
N17:=Nullspace(Matrix(AL17 - 1));

"Dimesion of eigenspace lambda=1 for w17 is:";
Dimension(N17);

N17c:=Nullspace(Matrix(AL17 + 1));

"Dimesion of eigenspace lambda=-1 for w17 is:";
Dimension(N17c);

B17:=[&+[(Integers()!(2*Eltseq(Basis(N17)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N17)]];
B17c:=[&+[(Integers()!(2*Eltseq(Basis(N17c)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N17c)]];

X119,jMap:= modformeqns(B17c cat B17, 119, 500, 17);
"Model for X0(119) is:";
X119;
"";
RR<[u]>:=CoordinateRing(AmbientSpace(X119));
n:=Dimension(AmbientSpace(X119));

H:=Matrix(RationalField(), 11, 11, [1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1]);
rows:=[[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]];
w17:=iso<X119->X119 | rows, rows>;
"w17 on X0(119) is given by:";
w17;
"";

X119w17, quotMap:=CurveQuotient(AutomorphismGroup(X119, [w17]));
"Genus of X0(119) is ", Genus(X119);
"Genus of X0(119)/w17 is ", Genus(X119w17);
"";

KK<rt>:=QuadraticField(-19);

P1:=X119![-1,0,0,0,0,0,0,1,0,0,0];
P2:=X119![1,0,0,0,0,0,0,1,0,0,0];
P3:=X119![-3/7,0,-2/7,-3/7,-2/7,-1/7,-1/7,1,0,0,0];
P4:=X119![3/7,0,2/7,3/7,2/7,1/7,1/7,1,0,0,0];

P5:=X119(KK)![1/7*(-2*rt + 1), 0, 1/7*(rt - 4), 1/14*(3*rt - 19), 1/7*(rt + 3), 1/14*(rt - 11), 1/14*(rt + 3), -2, 2, -1, 1];
P6:=X119(KK)![1/7*(-2*rt + 1), 0, 1/7*(rt - 4), 1/14*(3*rt - 19), 1/7*(rt + 3), 1/14*(rt - 11), 1/14*(rt + 3), 2, -2, 1, -1];

"By finding poles of j-map, we find that we have these 4 cusps:";
Poles(jMap);
"";
"and these 2 non-cusps (plus Galois conjugates):";
"P5: ", P5;
"P6: ", P6;
"";
"P5 and P6 are not pullbacks of rationals, they reduce via quotMap to:";
quotMap(P5), quotMap(P6);
"";

X3:=ChangeRing(X119,GF(3));
Cg,phi,psi:=ClassGroup(X3);
Z:=FreeAbelianGroup(1);
degr:=hom<Cg->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(Cg)]>;
JF3:=Kernel(degr); // This is isomorphic to J_X(\F_3).
"J0(119)(Q)_tor is isomorphic to a subgroup of:";
JF3;
"";

X5:=ChangeRing(X119,GF(5));
Cg,phi,psi:=ClassGroup(X5);
Z:=FreeAbelianGroup(1);
degr:=hom<Cg->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(Cg)]>;
JF5:=Kernel(degr); // This is isomorphic to J_X(\F_5).
"J0(119)(Q)_tor is isomorphic to a subgroup of:";
JF5;
"";

"By injecting torsion mod 3 and mod 5, we see that torsion has at most ", GCD(#JF3, #JF5), " elements.";
"";

"Clearly, cuspidal group is generated by Di := P1 - Pi, i in {2,3,4}";
"";

D2:=Divisor(P1)-Divisor(P2);
D3:=Divisor(P1)-Divisor(P3);
D4:=Divisor(P1)-Divisor(P4);

b, ff := IsPrincipal(288*D3);
"Is D3 of order 288?";
b and not(IsPrincipal(144*D3) or IsPrincipal(96*D3));
"";

Dt := 9*D4;
b, ff := IsPrincipal(8*Dt);
"Is Dt := 9*D4 of order 8?";
b and not(IsPrincipal(4*Dt));
"";
"Is 4*Dt of order 2 in J(Q)_tor/<D3>?";
not(IsPrincipal(4*Dt)) and not(IsPrincipal(4*Dt + 144*D3));

"Therefore, Dt := 9*D4 and D3 generate a torsion subgroup of order 288*8 = 2304, hence the whole torsion";

//known degree 1 places
pls1:=[Place(P1), Place(P2), Place(P3), Place(P4)];

//known degree 2 places, we know both of these are non pullbacks
pls2:=[Place(P5), Place(P6)];

deg2:=[1*Place(P5), 1*Place(P6)];
deg2pb:=[];

for i in [1..#pls1] do
	for j in [i..#pls1] do
		deg2:=Append(deg2, 1*pls1[i] + 1*pls1[j]);
		if w17(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb:=Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

deg2npb := [DD : DD in deg2 | not DD in deg2pb];

"We have found ", #deg2, " points on X_0(119)^2(Q).";
#deg2pb, "of them are pullbacks of rationals from X0(95)/w19.";
#deg2npb, "of them are non-pullbacks";

//Finally, we do the sieve.
A := AbelianGroup([8, 288]);
divs := [Dt, D3];
genusC := Genus(X119w17);
auts := [H];
bp := deg2pb[1];

//A == J(119)(Q)_tors, hence I = 1
I := 1;

primes := [5];

bret := MWSieve(deg2, primes, X119, A, divs, auts, genusC, deg2pb, deg2npb, I, bp);
"Succeeded in proving that we have all exceptional quadratic pts? (true if succeeded, number otherwise)";
bret;
assert bret eq true;

"This means we have found all exceptional quadratic points on X0(119).";
"Now we only need to find rational points on X0(119)/w17 and check their pullbacks.";

Hyp, mp := SimplifiedModel(X119w17);
"X0(119)/w17 is actually ", Hyp;

"Rank of J(X0(119)/w17)(Q) is", RankBound(Jacobian(Hyp));
"Hence rk(J(X0(119)/w17)(Q)) = rk(J0(119)(Q)) = 0.";
"";

ptsMax := GCD(#Jacobian(ChangeRing(Hyp, GF(3))), #Jacobian(ChangeRing(Hyp, GF(5))));
"J(X0(95)/w19)(Q)_tors has at most ", ptsMax, " points.";

Inf1 := PointsAtInfinity(Hyp)[1];
Inf2 := PointsAtInfinity(Hyp)[2];

Div1 := Inf2 - Inf1;

S := {};

for i := 1 to 9 do
	S := S join {i*Div1};
end for;

"We have found ", #S, " points on J(X0(119)/w17)(Q) out of possible ", ptsMax;
assert #S eq ptsMax;
"Hence, we know whole J(X0(119)/w17)(Q).";
"";

R:={};
for a in S do
	R:=R join SequenceToSet(Roots(a[1]));
end for;

"Using Mumford representations, we get that there are no rational points on X0(119)/w17";
"(except pts at infinity).";
"";

assert #R eq 0;

"Hence, only rational points are points at infinity.";
"";
for i in [1..#PointsAtInfinity(Hyp)] do
	"Pullback of point ", PointsAtInfinity(Hyp)[i], " is:";
	Decomposition(Pullback(quotMap*mp, Place(PointsAtInfinity(Hyp)[i])));
end for;

"Hence, the only quadratic point up to Galois conjugacy, apart from rational cusps, are:";
P5, P6;