//B. Vukorepa: This is part of code written by Ozman and Siksek and used by Box in https://arxiv.org/pdf/1906.05206.pdf.

// X is a projective curve over rationals,
// p prime of good reduction,
// D divisor on X,
// This reduces to a divisor on X/F_p.

reduce:=function(X,Xp,D);
	if Type(D) eq DivCrvElt then
		decomp:=Decomposition(D);
		return &+[ pr[2]*$$(X,Xp,pr[1]) : pr in decomp]; // Reduce the problem to reducing places.
	end if;
	R<[x]>:=CoordinateRing(AmbientSpace(X));
        assert Type(D) eq PlcCrvElt;
        if  (Degree(D) eq 1) and (#{Degree(xx) : xx in x} eq 1) then
		P:=D;
		m:=Rank(R);
		KX:=FunctionField(X);
		inds:=[i : i in [1..m] | &and[Valuation(KX!(x[j]/x[i]),P) ge 0 : j in [1..m]]];	
		assert #inds ne 0;
		i:=inds[1];
		PP:=[Evaluate(KX!(x[j]/x[i]),P) : j in [1..m]];
		denom:=LCM([Denominator(d) : d in PP]);
		PP:=[Integers()!(denom*d) : d in PP];
		g:=GCD(PP);
		PP:=[d div g : d in PP];
		Fp:=BaseRing(Xp);
		PP:=Xp![Fp!d : d in PP];
		return Place(PP);	
	end if;
	I:=Ideal(D);
	Fp:=BaseRing(Xp);
	p:=Characteristic(Fp);
	B:=Basis(I) cat DefiningEquations(X);
	m:=Rank(CoordinateRing(X));
	assert Rank(CoordinateRing(Xp)) eq m;
	R:=PolynomialRing(Integers(),m);
	BR:=[];
	for f in B do
		g:=f*p^-(Minimum([Valuation(c,p) : c in Coefficients(f)]));
		g:=g*LCM([Denominator(c) : c in Coefficients(g)]);
		Append(~BR,g);
	end for;
	J:=ideal<R | BR>;
	J:=Saturation(J,R!p);
	BR:=Basis(J);
	Rp:=CoordinateRing(AmbientSpace(Xp));
	assert Rank(Rp) eq m;
	BRp:=[Evaluate(f,[Rp.i : i in [1..m]]) : f in BR];
	Jp:=ideal<Rp| BRp>;
	Dp:=Divisor(Xp,Jp);
	return Dp;
end function;


