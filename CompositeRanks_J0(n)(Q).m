Ns := [60, 62, 69, 79, 83, 89, 92, 94, 95, 101, 119, 131];

for i in [1..#Ns] do
    "---------------";
    "Now determining r(J0(", Ns[i], ")(Q))...";
    M := CuspidalSubspace(ModularSymbols(Ns[i]));
    l := NewformDecomposition(M);

    rk := 0;

    "There are ", #l, " Galois orbits of Hecke eigenforms.";
    "";
    for j in [1..#l] do
        "Working with representative f_",j, "...";
        lRatio := LRatio(l[j], 1);
        if lRatio ne 0 then
            "LRatio is nonzero, hence, by Kolyvagin-Logachev r(A_f_",j,") = 0.";
        else 
            coef, ord := LSeriesLeadingCoefficient(l[j], 1, 20);
            "Order of vanishing of L(f_",j,", 1) is ", ord;
            "By Kolyvagin-Logachev: r(A_f_",j,") = ", Degree(HeckeEigenvalueField(l[j]));
            rk := rk + Degree(HeckeEigenvalueField(l[j]));
        end if;
        "";
  
    end for;
    "That proves that r(J0(", Ns[i], ")(Q)) = ",rk,".";
    "---------------";
    "";
end for;
