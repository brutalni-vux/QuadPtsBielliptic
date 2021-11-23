Magma V2.26-3     Sun Nov 21 2021 22:54:58 on euler    [Seed = 2021906812]
Type ? for help.  Type <Ctrl>-D to quit.
Loading "X0p_NiceModel.m"
Loading "Chabauty_MWSieve.m"
Loading "ChabautyHelp.m"
Dimension of CuspForms(79) is:  6
Dimesion of eigenspace lambda = 1 for w79 is:  1
Dimesion of eigenspace lambda = -1 for w79 is:  5

Nice model for X0(79) is:
Curve over Rational Field defined by
x[1]^2 - 2*x[1]*x[2] - x[2]^2 - 3*x[3]^2 - 2*x[3]*x[4] + 4*x[3]*x[5] + 3*x[4]^2 
    + 2*x[4]*x[5] - 11*x[5]^2 - x[6]^2,
x[1]*x[3] - x[2]^2 + 3*x[3]*x[4] - x[3]*x[5] - 2*x[4]^2 - x[4]*x[5] + 4*x[5]^2,
x[1]*x[4] - x[2]*x[3] + 2*x[3]*x[4] - x[4]^2 + x[5]^2,
x[1]*x[5] - x[3]^2 + 2*x[3]*x[4] + x[3]*x[5] - x[4]^2 - x[4]*x[5] + 2*x[5]^2,
x[2]*x[4] - x[3]^2 + 2*x[3]*x[5] - x[5]^2,
x[2]*x[5] - x[3]*x[4] + x[3]*x[5] + x[4]^2 - 2*x[5]^2

w79 on X0(79) is given by:
Mapping from: Crv: X79 to Crv: X79
with equations : 
u[1]
u[2]
u[3]
u[4]
u[5]
-u[6]
and inverse
u[1]
u[2]
u[3]
u[4]
u[5]
-u[6]

Genus of X0(79) is  6
Genus of X0(79)/w79 is  1

We have found these points on X0(79):
[ (1 : 0 : 0 : 0 : 0 : 1), (-1 : 0 : 0 : 0 : 0 : 1) ]

Is Dtor := pts[1] - pts[2] a generator for J0(79)(Q)_tors?  true

X0(79)/w79 is actually the following elliptic curve:
Elliptic Curve defined by y^2 - 25*x*y - 2*y = x^3 - 152*x^2 - 21*x over 
Rational Field

It has the following MW group:
Abelian Group isomorphic to Z
Defined on 1 generator (free)

We have found  3  points on X_0(79)^2(Q).
1 of them are pullbacks of rationals from X0(79)/w79.
3
The number of possible cosets unknown points can land in is 8
5
The number of possible cosets unknown points can land in is 0
MWSieve achieved its goal? (true if succeeded, number otherwise)
true
It follows that if there is an unknown Q in X0(79)^2(Q), then 2[Q - bp] is fixed
by w79.
This implies that [Q - bp] is fixed by w79 (since there is no 2-torsion)
Then Q ~ w79(Q), which implies that Q = w79(Q) because X0(79) is not 
hyperelliptic.
Then either Q is a pullback, or it is fixed by w79 pointwise.
If P = (X1:X2:X3:X4:X5:X6) is fixed by w79, then either X6 = 0 or P = 
(0:0:0:0:0:1) (not on X0(79))

We find all such P by imposing condition X6 = 0 and finding points on the 
scheme:

Scheme over Rational Field defined by
x[6],
x[1]^2 - 2*x[1]*x[2] - x[2]^2 - 3*x[3]^2 - 2*x[3]*x[4] + 4*x[3]*x[5] + 3*x[4]^2 
    + 2*x[4]*x[5] - 11*x[5]^2 - x[6]^2,
x[1]*x[3] - x[2]^2 + 3*x[3]*x[4] - x[3]*x[5] - 2*x[4]^2 - x[4]*x[5] + 4*x[5]^2,
x[1]*x[4] - x[2]*x[3] + 2*x[3]*x[4] - x[4]^2 + x[5]^2,
x[1]*x[5] - x[3]^2 + 2*x[3]*x[4] + x[3]*x[5] - x[4]^2 - x[4]*x[5] + 2*x[5]^2,
x[2]*x[4] - x[3]^2 + 2*x[3]*x[5] - x[5]^2,
x[2]*x[5] - x[3]*x[4] + x[3]*x[5] + x[4]^2 - 2*x[5]^2

All such P are:
{@ (-2/87*r1^4 - 37/87*r1^3 - 10/29*r1^2 - 13/87*r1 + 74/29 : -7/87*r1^4 + 
    1/87*r1^3 - 6/29*r1^2 + 85/87*r1 - 31/29 : 17/87*r1^4 + 10/87*r1^3 - 
    2/29*r1^2 - 107/87*r1 + 38/29 : r1 : 1 : 0), (-2/87*r2^4 - 37/87*r2^3 - 
    10/29*r2^2 - 13/87*r2 + 74/29 : -7/87*r2^4 + 1/87*r2^3 - 6/29*r2^2 + 
    85/87*r2 - 31/29 : 17/87*r2^4 + 10/87*r2^3 - 2/29*r2^2 - 107/87*r2 + 38/29 :
r2 : 1 : 0), (-2/87*r3^4 - 37/87*r3^3 - 10/29*r3^2 - 13/87*r3 + 74/29 : 
    -7/87*r3^4 + 1/87*r3^3 - 6/29*r3^2 + 85/87*r3 - 31/29 : 17/87*r3^4 + 
    10/87*r3^3 - 2/29*r3^2 - 107/87*r3 + 38/29 : r3 : 1 : 0), (-2/87*r4^4 - 
    37/87*r4^3 - 10/29*r4^2 - 13/87*r4 + 74/29 : -7/87*r4^4 + 1/87*r4^3 - 
    6/29*r4^2 + 85/87*r4 - 31/29 : 17/87*r4^4 + 10/87*r4^3 - 2/29*r4^2 - 
    107/87*r4 + 38/29 : r4 : 1 : 0), (-2/87*r5^4 - 37/87*r5^3 - 10/29*r5^2 - 
    13/87*r5 + 74/29 : -7/87*r5^4 + 1/87*r5^3 - 6/29*r5^2 + 85/87*r5 - 31/29 : 
    17/87*r5^4 + 10/87*r5^3 - 2/29*r5^2 - 107/87*r5 + 38/29 : r5 : 1 : 0), 
(2*r6^4 - 41/3*r6^3 + 62/3*r6^2 - 13/3*r6 - 2 : -r6^4 + 25/3*r6^3 - 46/3*r6^2 + 
    17/3*r6 + 1 : -r6^4 + 22/3*r6^3 - 10*r6^2 + 1/3*r6 + 2/3 : r6 : 1 : 0), 
(2*r7^4 - 41/3*r7^3 + 62/3*r7^2 - 13/3*r7 - 2 : -r7^4 + 25/3*r7^3 - 46/3*r7^2 + 
    17/3*r7 + 1 : -r7^4 + 22/3*r7^3 - 10*r7^2 + 1/3*r7 + 2/3 : r7 : 1 : 0), 
(2*r8^4 - 41/3*r8^3 + 62/3*r8^2 - 13/3*r8 - 2 : -r8^4 + 25/3*r8^3 - 46/3*r8^2 + 
    17/3*r8 + 1 : -r8^4 + 22/3*r8^3 - 10*r8^2 + 1/3*r8 + 2/3 : r8 : 1 : 0), 
(2*r9^4 - 41/3*r9^3 + 62/3*r9^2 - 13/3*r9 - 2 : -r9^4 + 25/3*r9^3 - 46/3*r9^2 + 
    17/3*r9 + 1 : -r9^4 + 22/3*r9^3 - 10*r9^2 + 1/3*r9 + 2/3 : r9 : 1 : 0), 
(2*r10^4 - 41/3*r10^3 + 62/3*r10^2 - 13/3*r10 - 2 : -r10^4 + 25/3*r10^3 - 
    46/3*r10^2 + 17/3*r10 + 1 : -r10^4 + 22/3*r10^3 - 10*r10^2 + 1/3*r10 + 2/3 :
r10 : 1 : 0) @}

Clearly, none such points are quadratic.
Hence there are no quadratic points on X0(79) not coming from pullbacks of 
rationals.

Total time: 43.079 seconds, Total memory usage: 184.72MB
