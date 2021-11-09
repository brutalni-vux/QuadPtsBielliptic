//see https://wstein.org/edu/Fall2003/252/references/magma/ModSym.pdf
//for Modular Symbols package

Composite_Ns := [60,62,69,92,94,95,119];

for i in [1..#Composite_Ns] do
    "---------------";
    "Now proving that r(J0(", Composite_Ns[i], ")(Q)) is 0...";
    M := CuspidalSubspace(ModularSymbols(Composite_Ns[i]));
    l:=NewformDecomposition(M);

    "There are ", #l, " Galois orbits of Hecke eigenforms.";
    "Values of L(A_f, 1)/Omega_{A_f} for representatives f are:";
    for j in [1..#l] do
        LRatio(l[j], 1);
    end for;
    "All values are nonzero!";
    "Hence, by Kolyvagin-Logachev r(A_f) = 0 for all f";
    "That proves that r(J0(", Composite_Ns[i], ")(Q) = 0.";
    "---------------";
    "";
end for;