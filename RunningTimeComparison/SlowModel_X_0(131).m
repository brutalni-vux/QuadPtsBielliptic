load "X0p_NiceModel.m";

//we compute the model for X_0(131) using the default cuspform basis in Magma

C := CuspForms(131);
"Dimension of CuspForms(131) is: ", Dimension(C);

X131 := modformeqns(Basis(C), 131, 500, 1);
"Nice model for X0(131) is:";
X131;
