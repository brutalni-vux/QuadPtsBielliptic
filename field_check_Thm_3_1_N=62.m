//The code in this file checks the field of definition of the CM points on X0(62)

//we use the version used for prime levels X0(p) because it has Curve_and_Map function
//there is no other difference and the model coputations are the same as with non-primes
//we also don't need the jMap
load "X0p_NiceModel.m";

our_eq := function(n, d)
	C := CuspForms(n);
	AL := AtkinLehnerOperator(C, d);

	N := Nullspace(Matrix(AL- 1));
	Nc := Nullspace(Matrix(AL + 1));

	B := [&+[(Integers()!(2*Eltseq(Basis(N)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N)]];
	Bc := [&+[(Integers()!(2*Eltseq(Basis(Nc)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(Nc)]];

	X := modformeqns(Bc cat B, n, 500, 1);

	Xwd, mp:=Curve_and_Map(X, Dimension(Nc));
	assert Genus(Xwd) eq 1;

	pts := PointSearch(Xwd, 2);
	E, f := EllipticCurve(Xwd, pts[1]);
	return X, E, mp*f;
end function;


///// check for 62: ////
X, Xwd, mp := our_eq(62, 31);
assert Rank(Xwd) eq 0;

g, m := TorsionSubgroup(Xwd);
P := m(g.1);

K := ResidueClassField(Decomposition(Pullback(mp, Place(P)))[1, 1]);
assert SquareFreeFactorization(Discriminant(Integers(K))) eq -3;

K := ResidueClassField(Decomposition(Pullback(mp,Place(4*P)))[1, 1]);
assert SquareFreeFactorization(Discriminant(Integers(K))) eq -3;
