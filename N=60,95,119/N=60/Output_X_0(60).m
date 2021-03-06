Magma V2.26-3     Fri Nov 12 2021 18:27:36 on euler    [Seed = 903540855]
Type ? for help.  Type <Ctrl>-D to quit.
Loading "quadptssieve.m"
Loading "ozmansiksek.m"
Loading "ozmansiksek_TorsHelp.m"
Loading "X0N_NiceModel.m"
Loading "HelpTorsion.m"
Dimension of CuspForms(60) is:
7
Dimesion of eigenspace lambda=1 for w15 is:
1
Dimesion of eigenspace lambda=-1 for w15 is:
6
Model for X0(60) is:
Curve over Rational Field defined by
x[1]^2 + 2*x[2]^2 - x[3]^2 + 6*x[4]^2 - 6*x[4]*x[6] + x[5]^2 + 4*x[6]^2 - 
    x[7]^2,
x[1]*x[3] - x[2]^2 + x[4]^2 + x[4]*x[6],
x[1]*x[4] - x[2]*x[3],
x[1]*x[5] - x[3]^2 - x[4]^2 - x[4]*x[6] + x[6]^2,
x[1]*x[6] - x[3]*x[4] - x[4]*x[5],
x[2]*x[4] - x[3]^2 - x[4]*x[6],
x[2]*x[5] - x[3]*x[4] - x[4]*x[5],
x[2]*x[6] - x[4]^2 - x[4]*x[6],
x[3]*x[5] - x[4]^2 - x[4]*x[6] + x[6]^2,
x[3]*x[6] - x[4]*x[5]

w15 on X0(60) is given by:
Mapping from: Crv: X60 to Crv: X60
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

Genus of X0(60) is  7
Genus of X0(60)/w15 is  1

By finding poles of j-map, we find that we have these 12 cusps (P1, ..., P12):
[
    Place at (0 : 0 : 1/4 : 1/4 : -1/4 : -1/4 : 1),
    Place at (-1/2 : -1/2 : -1/4 : -1/4 : -1/4 : -1/4 : 1),
    Place at (0 : 0 : -1/4 : -1/4 : 1/4 : 1/4 : 1),
    Place at (1/2 : 1/2 : 1/4 : 1/4 : 1/4 : 1/4 : 1),
    Place at (0 : 0 : -1/4 : 1/4 : 1/4 : -1/4 : 1),
    Place at (1/2 : -1/2 : 1/4 : -1/4 : 1/4 : -1/4 : 1),
    Place at (0 : 0 : 1/4 : -1/4 : -1/4 : 1/4 : 1),
    Place at (-1/2 : 1/2 : -1/4 : 1/4 : -1/4 : 1/4 : 1),
    Place at (-1 : 0 : 0 : 0 : 0 : 0 : 1),
    Place at (1 : 0 : 0 : 0 : 0 : 0 : 1),
    Place at (0 : 0 : 0 : 0 : 1 : 0 : 1),
    Place at (0 : 0 : 0 : 0 : -1 : 0 : 1)
]

Cuspidal group is generated by Di := P1 - Pi, i in {2,3, ..., 12}

Is D2 of order 24?
true

Is D3 of order 24?
true

Is a*D2 = b*D3 while (a,b)!=(0,0) mod 24? (are D2 and D3 linearly dependent?)
false , hence, D2 and D3 generate subgroup iso. to (Z/24Z)^2

Is D5 of order 24?
true

Is coset of D5 of order 24 in J(60)(Q)_tor/<D2, D3>?
true , hence, we know that D2, D3, D5 generate a subgroup (Z/24Z)^3

Is Dt := D7 - 4*D2 - 6*D3 of order 4?
true

Is coset of Dt of order 4 in J(60)(Q)_tor/<D2,D3,D5>?
true
Hence we have proven that Dt, D5, D3, D2 generate a subgroup Z/4Z + (Z/24Z)^3

Now we prove that J0(60)(Q)_tors == Z/4Z + (Z/24Z)^3

By Ozman, Siksek, cuspidal subgroup is isomorphic to:
Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/24
Defined on 4 generators in supergroup:
    Ksub.1 = $.3
    Ksub.2 = $.4
    Ksub.3 = 2*$.5
    Ksub.4 = 6*$.6
Relations:
    4*Ksub.1 = 0
    24*Ksub.2 = 0
    24*Ksub.3 = 0
    24*Ksub.4 = 0
(same as we got from our computations)

J0(60)(Q)_tors is isomorphic to one of the following groups:

{
    Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/24
    Defined on 4 generators in supergroup:
        $.1 = $.2 + $.3
        $.2 = $.4
        $.3 = 2*$.5
        $.4 = 6*$.6
    Relations:
        4*$.1 = 0
        24*$.2 = 0
        24*$.3 = 0
        24*$.4 = 0,

    Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/48
    Defined on 4 generators in supergroup:
        $.1 = $.2 + $.3
        $.2 = $.4
        $.3 = 2*$.5
        $.4 = $.2 + 3*$.6
    Relations:
        4*$.1 = 0
        24*$.2 = 0
        24*$.3 = 0
        48*$.4 = 0,

    Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/48
    Defined on 4 generators in supergroup:
        $.1 = $.2 + $.3
        $.2 = $.4
        $.3 = 2*$.5
        $.4 = $.1 + $.3 + 3*$.6
    Relations:
        4*$.1 = 0
        24*$.2 = 0
        24*$.3 = 0
        48*$.4 = 0,

    Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/48
    Defined on 4 generators in supergroup:
        $.1 = $.2 + $.3
        $.2 = $.4
        $.3 = 6*$.6
        $.4 = $.2 + $.5 + 3*$.6
    Relations:
        4*$.1 = 0
        24*$.2 = 0
        24*$.3 = 0
        48*$.4 = 0,

    Abelian Group isomorphic to Z/4 + Z/24 + Z/24 + Z/48
    Defined on 4 generators in supergroup:
        $.1 = $.2 + $.3
        $.2 = $.4
        $.3 = 6*$.6
        $.4 = $.1 + $.3 + $.5 + 3*$.6
    Relations:
        4*$.1 = 0
        24*$.2 = 0
        24*$.3 = 0
        48*$.4 = 0
}

J(F23) is iso to:
Abelian Group isomorphic to Z/12 + Z/24 + Z/24 + Z/24 + Z/72 + Z/216
Defined on 6 generators in supergroup CG:
    JF23.1 = CG.1
    JF23.2 = CG.2
    JF23.3 = CG.3
    JF23.4 = CG.4
    JF23.5 = CG.5
    JF23.6 = CG.6
Relations:
    12*JF23.1 = 0
    24*JF23.2 = 0
    24*JF23.3 = 0
    24*JF23.4 = 0
    72*JF23.5 = 0
    216*JF23.6 = 0

This shows that J0(60)(Q)_tors can't have an element of order 48.
By the previous information and with rk(J0(60)(Q)) = 0, we conclude J0(60)(Q) ==
Z/4Z + (Z/24Z)^3

We have found  78  points on X_0(60)^2(Q).
6 of them are pullbacks of rationals from X0(60)/w15.
72 of them are non-pullbacks
13
13
Succeeded in proving that we have found all exceptional quadratic pts? (true if 
succeeded, number otherwise)
true
Hence, there are no quadratic points on X0(60) not coming from X0(60)/w15(Q).
We now have to find rational points on X0(60)/w15(Q) and check their pullbacks.

X0(60)/w15(Q) is actually the following elliptic curve E/Q:
Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field

It has rank 0
Notice that this shows that rk(J0(60)(Q)) == rk(J(C)(Q)).

Its torsion is:
Abelian Group isomorphic to Z/6
Defined on 1 generator
Relations:
    6*grp.1 = 0

Here are all rational points on E:
{@ (0 : 1 : 0), (-1 : 1 : 1), (-1 : -1 : 1), (0 : 0 : 1), (1 : 1 : 1), (1 : -1 :
1) @}

Pullback of point  (0 : 1 : 0)  is:
[
    <Place at (0 : 0 : 0 : 0 : 1 : 0 : 1), 1>,
    <Place at (0 : 0 : 0 : 0 : -1 : 0 : 1), 1>
]

Pullback of point  (-1 : 1 : 1)  is:
[
    <Place at (-1/2 : -1/2 : -1/4 : -1/4 : -1/4 : -1/4 : 1), 1>,
    <Place at (1/2 : 1/2 : 1/4 : 1/4 : 1/4 : 1/4 : 1), 1>
]

Pullback of point  (-1 : -1 : 1)  is:
[
    <Place at (1/2 : -1/2 : 1/4 : -1/4 : 1/4 : -1/4 : 1), 1>,
    <Place at (-1/2 : 1/2 : -1/4 : 1/4 : -1/4 : 1/4 : 1), 1>
]

Pullback of point  (0 : 0 : 1)  is:
[
    <Place at (-1 : 0 : 0 : 0 : 0 : 0 : 1), 1>,
    <Place at (1 : 0 : 0 : 0 : 0 : 0 : 1), 1>
]

Pullback of point  (1 : 1 : 1)  is:
[
    <Place at (0 : 0 : -1/4 : 1/4 : 1/4 : -1/4 : 1), 1>,
    <Place at (0 : 0 : 1/4 : -1/4 : -1/4 : 1/4 : 1), 1>
]

Pullback of point  (1 : -1 : 1)  is:
[
    <Place at (0 : 0 : 1/4 : 1/4 : -1/4 : -1/4 : 1), 1>,
    <Place at (0 : 0 : -1/4 : -1/4 : 1/4 : 1/4 : 1), 1>
]

Hence, there are no new quadratic points on X0(60) coming from X0(60)/w15(Q).

Total time: 1395.940 seconds, Total memory usage: 1736.44MB
