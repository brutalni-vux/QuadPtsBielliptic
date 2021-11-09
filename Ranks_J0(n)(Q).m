//this code was written by us, although some of this info is already known
//see https://wstein.org/edu/Fall2003/252/references/magma/ModSym.pdf for info about modular symbols package
//if possible, this code determines rk(J0(N)(Q)) for all N in Ns using Kolyvagin-Logachev theorem

Ns := [60, 62, 69, 79, 83, 89, 92, 94, 95, 101, 119, 131];

for i in [1..#Ns] do
    "---------------";
    "Now determining r(J0(", Ns[i], ")(Q))...";
    M := CuspidalSubspace(ModularSymbols(Ns[i]));
    l := NewformDecomposition(M);

    rk := 0;
    provable := true;

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
            if ord eq 1 then
                "By Kolyvagin-Logachev: r(A_f_",j,") = ", Degree(HeckeEigenvalueField(l[j]));
                rk := rk + Degree(HeckeEigenvalueField(l[j]));
            else
                "We can't use Kolyvagin-Logachev to determine r(A_f_",j,").";
                provable := false;
            end if;
        end if;
        "";
  
    end for;
    if provable then
        "That proves that r(J0(", Ns[i], ")(Q)) = ",rk,".";
    else
        "We were unable to provably determine r(J0(", Ns[i], ")(Q))";
    end if;
    "---------------";
    "";
end for;
