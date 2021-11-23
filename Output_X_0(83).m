Magma V2.26-3     Sun Nov 21 2021 22:54:44 on euler    [Seed = 434288651]
Type ? for help.  Type <Ctrl>-D to quit.
Loading "X0p_NiceModel.m"
Loading "Chabauty_MWSieve.m"
Loading "ChabautyHelp.m"
Dimension of CuspForms(83) is:  7
Dimesion of eigenspace lambda = 1 for w83 is:  1
Dimesion of eigenspace lambda = -1 for w83 is:  6

Nice model for X0(83) is:
Curve over Rational Field defined by
x[1]^2 - 2*x[1]*x[2] - x[2]^2 - x[3]^2 + 2*x[3]*x[4] - x[4]^2 - 4*x[4]*x[5] + 
    2*x[4]*x[6] - 10*x[5]^2 + 24*x[5]*x[6] - 31*x[6]^2 - x[7]^2,
x[1]*x[3] - x[2]^2 + 2*x[3]*x[4] - 2*x[4]^2 + 6*x[4]*x[5] - 6*x[4]*x[6] - 
    3*x[5]^2 + 6*x[5]*x[6],
x[1]*x[4] - x[2]*x[3] + 2*x[3]*x[4] - 2*x[4]^2 + 5*x[4]*x[5] - 2*x[4]*x[6] - 
    4*x[5]^2 + 8*x[5]*x[6] - 6*x[6]^2,
x[1]*x[5] - x[3]^2 + x[3]*x[4] + x[4]*x[5] - 3*x[4]*x[6] + 4*x[5]^2 - 
    6*x[5]*x[6] + 9*x[6]^2,
x[1]*x[6] - x[4]^2 + 2*x[4]*x[5] - x[5]^2 + 2*x[5]*x[6],
x[2]*x[4] - x[3]^2 + 2*x[4]^2 - 4*x[4]*x[5] - x[4]*x[6] + 7*x[5]^2 - 
    14*x[5]*x[6] + 15*x[6]^2,
x[2]*x[5] - x[3]*x[4] + x[4]^2 - 2*x[4]*x[5] + 2*x[4]*x[6] + 3*x[5]^2 - 
    8*x[5]*x[6] + 6*x[6]^2,
x[2]*x[6] - x[4]*x[5] + x[4]*x[6] + x[5]^2 - 2*x[5]*x[6],
x[3]*x[5] - x[4]^2 + x[4]*x[5] + x[4]*x[6] - 3*x[5]^2 + 5*x[5]*x[6] - 6*x[6]^2,
x[3]*x[6] - x[5]^2 + 2*x[5]*x[6] - 3*x[6]^2

w83 on X0(83) is given by:
Mapping from: Crv: X83 to Crv: X83
with equations : 
u[1]
u[2]
u[3]
u[4]
u[5]
u[6]
-u[7]
and inverse
u[1]
u[2]
u[3]
u[4]
u[5]
u[6]
-u[7]

Genus of X0(83) is  7
Genus of X0(83)/w83 is  1

We have found these points on X0(83):
[ (1 : 0 : 0 : 0 : 0 : 0 : 1), (-1 : 0 : 0 : 0 : 0 : 0 : 1) ]

Is Dtor := pts[1] - pts[2] a generator for J0(83)(Q)_tors?  true

X0(83)/w83 is actually the following elliptic curve:
Elliptic Curve defined by y^2 + 92*x*y - 59904*y = x^3 - 7556*x^2 + 12684288*x -
6979911680 over Rational Field

It has the following MW group:
Abelian Group isomorphic to Z
Defined on 1 generator (free)

We have found  3  points on X_0(83)^2(Q).
1 of them are pullbacks of rationals from X0(83)/w83.
3
The number of possible cosets unknown points can land in is 2
5
The number of possible cosets unknown points can land in is 0
MWSieve achieved its goal? (true if succeeded, number otherwise)
true
It follows that if there is an unknown Q in X0(83)^2(Q), then 2[Q - bp] is fixed
by w83.
This implies that [Q - bp] is fixed by w83 (since there is no 2-torsion)
Then Q ~ w83(Q), which implies that Q = w83(Q) because X0(83) is not 
hyperelliptic.
Then either Q is a pullback, or it is fixed by w83 pointwise.
If P = (X1:X2:X3:X4:X5:X6:X7) is fixed by w83, then either X7 = 0 or P = 
(0:0:0:0:0:0:1) (not on X0(83))

We find all such P by imposing condition X7 = 0 and finding points on the 
scheme:

Scheme over Rational Field defined by
x[7],
x[1]^2 - 2*x[1]*x[2] - x[2]^2 - x[3]^2 + 2*x[3]*x[4] - x[4]^2 - 4*x[4]*x[5] + 
    2*x[4]*x[6] - 10*x[5]^2 + 24*x[5]*x[6] - 31*x[6]^2 - x[7]^2,
x[1]*x[3] - x[2]^2 + 2*x[3]*x[4] - 2*x[4]^2 + 6*x[4]*x[5] - 6*x[4]*x[6] - 
    3*x[5]^2 + 6*x[5]*x[6],
x[1]*x[4] - x[2]*x[3] + 2*x[3]*x[4] - 2*x[4]^2 + 5*x[4]*x[5] - 2*x[4]*x[6] - 
    4*x[5]^2 + 8*x[5]*x[6] - 6*x[6]^2,
x[1]*x[5] - x[3]^2 + x[3]*x[4] + x[4]*x[5] - 3*x[4]*x[6] + 4*x[5]^2 - 
    6*x[5]*x[6] + 9*x[6]^2,
x[1]*x[6] - x[4]^2 + 2*x[4]*x[5] - x[5]^2 + 2*x[5]*x[6],
x[2]*x[4] - x[3]^2 + 2*x[4]^2 - 4*x[4]*x[5] - x[4]*x[6] + 7*x[5]^2 - 
    14*x[5]*x[6] + 15*x[6]^2,
x[2]*x[5] - x[3]*x[4] + x[4]^2 - 2*x[4]*x[5] + 2*x[4]*x[6] + 3*x[5]^2 - 
    8*x[5]*x[6] + 6*x[6]^2,
x[2]*x[6] - x[4]*x[5] + x[4]*x[6] + x[5]^2 - 2*x[5]*x[6],
x[3]*x[5] - x[4]^2 + x[4]*x[5] + x[4]*x[6] - 3*x[5]^2 + 5*x[5]*x[6] - 6*x[6]^2,
x[3]*x[6] - x[5]^2 + 2*x[5]*x[6] - 3*x[6]^2

All such P are:
{@ (118/535*r1^8 - 2518/535*r1^7 + 17401/535*r1^6 - 61862/535*r1^5 + 
    26096/107*r1^4 - 172294/535*r1^3 + 142047/535*r1^2 - 68182/535*r1 + 2668/107
: -118/535*r1^8 + 2518/535*r1^7 - 17401/535*r1^6 + 61862/535*r1^5 - 
    26096/107*r1^4 + 172829/535*r1^3 - 144722/535*r1^2 + 72462/535*r1 - 3310/107
: r1^2 - 2*r1 + 3 : -27/535*r1^8 + 422/535*r1^7 - 894/535*r1^6 - 4552/535*r1^5 +
    5360/107*r1^4 - 58509/535*r1^3 + 68717/535*r1^2 - 44582/535*r1 + 2770/107 : 
r1 : 1 : 0), (118/535*r2^8 - 2518/535*r2^7 + 17401/535*r2^6 - 61862/535*r2^5 + 
    26096/107*r2^4 - 172294/535*r2^3 + 142047/535*r2^2 - 68182/535*r2 + 2668/107
: -118/535*r2^8 + 2518/535*r2^7 - 17401/535*r2^6 + 61862/535*r2^5 - 
    26096/107*r2^4 + 172829/535*r2^3 - 144722/535*r2^2 + 72462/535*r2 - 3310/107
: r2^2 - 2*r2 + 3 : -27/535*r2^8 + 422/535*r2^7 - 894/535*r2^6 - 4552/535*r2^5 +
    5360/107*r2^4 - 58509/535*r2^3 + 68717/535*r2^2 - 44582/535*r2 + 2770/107 : 
r2 : 1 : 0), (118/535*r3^8 - 2518/535*r3^7 + 17401/535*r3^6 - 61862/535*r3^5 + 
    26096/107*r3^4 - 172294/535*r3^3 + 142047/535*r3^2 - 68182/535*r3 + 2668/107
: -118/535*r3^8 + 2518/535*r3^7 - 17401/535*r3^6 + 61862/535*r3^5 - 
    26096/107*r3^4 + 172829/535*r3^3 - 144722/535*r3^2 + 72462/535*r3 - 3310/107
: r3^2 - 2*r3 + 3 : -27/535*r3^8 + 422/535*r3^7 - 894/535*r3^6 - 4552/535*r3^5 +
    5360/107*r3^4 - 58509/535*r3^3 + 68717/535*r3^2 - 44582/535*r3 + 2770/107 : 
r3 : 1 : 0), (118/535*r4^8 - 2518/535*r4^7 + 17401/535*r4^6 - 61862/535*r4^5 + 
    26096/107*r4^4 - 172294/535*r4^3 + 142047/535*r4^2 - 68182/535*r4 + 2668/107
: -118/535*r4^8 + 2518/535*r4^7 - 17401/535*r4^6 + 61862/535*r4^5 - 
    26096/107*r4^4 + 172829/535*r4^3 - 144722/535*r4^2 + 72462/535*r4 - 3310/107
: r4^2 - 2*r4 + 3 : -27/535*r4^8 + 422/535*r4^7 - 894/535*r4^6 - 4552/535*r4^5 +
    5360/107*r4^4 - 58509/535*r4^3 + 68717/535*r4^2 - 44582/535*r4 + 2770/107 : 
r4 : 1 : 0), (118/535*r5^8 - 2518/535*r5^7 + 17401/535*r5^6 - 61862/535*r5^5 + 
    26096/107*r5^4 - 172294/535*r5^3 + 142047/535*r5^2 - 68182/535*r5 + 2668/107
: -118/535*r5^8 + 2518/535*r5^7 - 17401/535*r5^6 + 61862/535*r5^5 - 
    26096/107*r5^4 + 172829/535*r5^3 - 144722/535*r5^2 + 72462/535*r5 - 3310/107
: r5^2 - 2*r5 + 3 : -27/535*r5^8 + 422/535*r5^7 - 894/535*r5^6 - 4552/535*r5^5 +
    5360/107*r5^4 - 58509/535*r5^3 + 68717/535*r5^2 - 44582/535*r5 + 2770/107 : 
r5 : 1 : 0), (118/535*r6^8 - 2518/535*r6^7 + 17401/535*r6^6 - 61862/535*r6^5 + 
    26096/107*r6^4 - 172294/535*r6^3 + 142047/535*r6^2 - 68182/535*r6 + 2668/107
: -118/535*r6^8 + 2518/535*r6^7 - 17401/535*r6^6 + 61862/535*r6^5 - 
    26096/107*r6^4 + 172829/535*r6^3 - 144722/535*r6^2 + 72462/535*r6 - 3310/107
: r6^2 - 2*r6 + 3 : -27/535*r6^8 + 422/535*r6^7 - 894/535*r6^6 - 4552/535*r6^5 +
    5360/107*r6^4 - 58509/535*r6^3 + 68717/535*r6^2 - 44582/535*r6 + 2770/107 : 
r6 : 1 : 0), (118/535*r7^8 - 2518/535*r7^7 + 17401/535*r7^6 - 61862/535*r7^5 + 
    26096/107*r7^4 - 172294/535*r7^3 + 142047/535*r7^2 - 68182/535*r7 + 2668/107
: -118/535*r7^8 + 2518/535*r7^7 - 17401/535*r7^6 + 61862/535*r7^5 - 
    26096/107*r7^4 + 172829/535*r7^3 - 144722/535*r7^2 + 72462/535*r7 - 3310/107
: r7^2 - 2*r7 + 3 : -27/535*r7^8 + 422/535*r7^7 - 894/535*r7^6 - 4552/535*r7^5 +
    5360/107*r7^4 - 58509/535*r7^3 + 68717/535*r7^2 - 44582/535*r7 + 2770/107 : 
r7 : 1 : 0), (118/535*r8^8 - 2518/535*r8^7 + 17401/535*r8^6 - 61862/535*r8^5 + 
    26096/107*r8^4 - 172294/535*r8^3 + 142047/535*r8^2 - 68182/535*r8 + 2668/107
: -118/535*r8^8 + 2518/535*r8^7 - 17401/535*r8^6 + 61862/535*r8^5 - 
    26096/107*r8^4 + 172829/535*r8^3 - 144722/535*r8^2 + 72462/535*r8 - 3310/107
: r8^2 - 2*r8 + 3 : -27/535*r8^8 + 422/535*r8^7 - 894/535*r8^6 - 4552/535*r8^5 +
    5360/107*r8^4 - 58509/535*r8^3 + 68717/535*r8^2 - 44582/535*r8 + 2770/107 : 
r8 : 1 : 0), (118/535*r9^8 - 2518/535*r9^7 + 17401/535*r9^6 - 61862/535*r9^5 + 
    26096/107*r9^4 - 172294/535*r9^3 + 142047/535*r9^2 - 68182/535*r9 + 2668/107
: -118/535*r9^8 + 2518/535*r9^7 - 17401/535*r9^6 + 61862/535*r9^5 - 
    26096/107*r9^4 + 172829/535*r9^3 - 144722/535*r9^2 + 72462/535*r9 - 3310/107
: r9^2 - 2*r9 + 3 : -27/535*r9^8 + 422/535*r9^7 - 894/535*r9^6 - 4552/535*r9^5 +
    5360/107*r9^4 - 58509/535*r9^3 + 68717/535*r9^2 - 44582/535*r9 + 2770/107 : 
r9 : 1 : 0), (4*r10^2 - 17*r10 + 20 : -2*r10^2 + 8*r10 - 10 : r10^2 - 2*r10 + 3 
: -r10^2 + 5*r10 - 6 : r10 : 1 : 0), (4*r11^2 - 17*r11 + 20 : -2*r11^2 + 8*r11 -
    10 : r11^2 - 2*r11 + 3 : -r11^2 + 5*r11 - 6 : r11 : 1 : 0), (4*r12^2 - 
    17*r12 + 20 : -2*r12^2 + 8*r12 - 10 : r12^2 - 2*r12 + 3 : -r12^2 + 5*r12 - 6
: r12 : 1 : 0) @}

Clearly, none such points are quadratic.
Hence there are no quadratic points on X0(83) not coming from pullbacks of 
rationals.

Total time: 53.130 seconds, Total memory usage: 1704.41MB
