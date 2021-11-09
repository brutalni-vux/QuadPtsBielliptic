//see https://wstein.org/edu/Fall2003/252/references/magma/ModSym.pdf
//for Modular Symbols package

Rank0_Ns := [60,62,69,94,95,119];
Rank1_N := 92;

for i in [1..#Rank0_Ns] do
    "---------------";
    "Now proving that r(J0(", Rank0_Ns[i], ")(Q)) is 0...";
    M := CuspidalSubspace(ModularSymbols(Rank0_Ns[i]));
    l:=NewformDecomposition(M);

    "There are ", #l, " Galois orbits of Hecke eigenforms.";
    "Values of L(A_f, 1)/Omega_{A_f} for representatives f are:";
    for j in [1..#l] do
        LRatio(l[j], 1);
    end for;
    "All values are nonzero!";
    "Hence, by Kolyvagin-Logachev r(A_f) = 0 for all f";
    "That proves that r(J0(", Rank0_Ns[i], ")(Q)) = 0.";
    "---------------";
    "";
end for;

"---------------";
"Now proving that r(J0(", Rank1_N, ")(Q)) is 1...";
M := CuspidalSubspace(ModularSymbols(Rank1_N));
l := NewformDecomposition(M);

"There are ", #l, " Galois orbits of Hecke eigenforms.";
"Values of L(A_f, 1)/Omega_{A_f} for representatives f are:";
for j in [1..#l] do
    LRatio(l[j], 1);
end for;
"Hence, by Kolyvagin-Logachev r(A_f_i) = 0 for all i != 2.";

coef, ord := LSeriesLeadingCoefficient(l[2], 1, 20);
"Order of vanishing of L(f_2, 1) is ", ord;

Degree(HeckeEigenvalueField(l[2]));

"Again by Kolyvagin-Logachev: r(A_f_2) = ", Degree(HeckeEigenvalueField(l[2]));
"That proves that r(J0(", Rank1_N, ")(Q)) = 1.";
"---------------";
"";
