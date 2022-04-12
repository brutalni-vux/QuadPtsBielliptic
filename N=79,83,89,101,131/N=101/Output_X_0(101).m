Magma V2.26-3     Tue Apr 12 2022 21:28:54 on euler    [Seed = 2561615932]
Type ? for help.  Type <Ctrl>-D to quit.
Loading "X0p_NiceModel.m"
Loading "Chabauty_MWSieve.m"
Loading "ChabautyHelp.m"
Dimension of CuspForms(101) is:  8
Dimesion of eigenspace lambda = 1 for w101 is:  1
Dimesion of eigenspace lambda = -1 for w101 is:  7

Nice model for X0(101) is:
Curve over Rational Field defined by
x[1]^2 - 4*x[2]^2 - 4*x[2]*x[3] + 2*x[3]^2 + 8*x[3]*x[4] + 4*x[4]^2 - 
    20*x[4]*x[5] - 41*x[5]^2 + 112*x[5]*x[6] + 18*x[5]*x[7] - 22*x[6]^2 - 
    204*x[6]*x[7] + 170*x[7]^2 - x[8]^2,
x[1]*x[3] - x[2]^2 + 8*x[4]*x[5] - 8*x[5]^2 - 15*x[5]*x[6] + 21*x[5]*x[7] + 
    32*x[6]^2 - 61*x[6]*x[7] + 31*x[7]^2,
x[1]*x[4] - x[2]*x[3] + 8*x[4]*x[5] - 10*x[5]^2 - 5*x[5]*x[6] + 17*x[5]*x[7] + 
    22*x[6]^2 - 54*x[6]*x[7] + 30*x[7]^2,
x[1]*x[5] - x[3]^2 + 6*x[4]*x[5] - x[5]^2 - 20*x[5]*x[6] + 15*x[5]*x[7] + 
    25*x[6]^2 - 30*x[6]*x[7] + 9*x[7]^2,
x[1]*x[6] - x[3]*x[4] + 4*x[4]*x[5] - x[5]^2 - 10*x[5]*x[6] + 6*x[5]*x[7] + 
    12*x[6]^2 - 8*x[6]*x[7],
x[1]*x[7] - x[4]^2 + 2*x[4]*x[5] - x[5]*x[6] - x[5]*x[7] - 2*x[6]^2 + 
    10*x[6]*x[7] - 5*x[7]^2,
x[2]*x[4] - x[3]^2 + 8*x[5]^2 - 18*x[5]*x[6] - 2*x[5]*x[7] + 5*x[6]^2 + 
    26*x[6]*x[7] - 23*x[7]^2,
x[2]*x[5] - x[3]*x[4] + 6*x[5]^2 - 9*x[5]*x[6] - 6*x[5]*x[7] - 4*x[6]^2 + 
    33*x[6]*x[7] - 24*x[7]^2,
x[2]*x[6] - x[4]^2 + 4*x[5]^2 - x[5]*x[6] - 8*x[5]*x[7] - 10*x[6]^2 + 
    32*x[6]*x[7] - 20*x[7]^2,
x[2]*x[7] - x[4]*x[5] + 2*x[5]^2 - 3*x[5]*x[7] - 4*x[6]^2 + 12*x[6]*x[7] - 
    8*x[7]^2,
x[3]*x[5] - x[4]^2 + 6*x[5]*x[6] - 7*x[5]*x[7] - 10*x[6]^2 + 18*x[6]*x[7] - 
    8*x[7]^2,
x[3]*x[6] - x[4]*x[5] + 4*x[5]*x[6] - 2*x[5]*x[7] - 5*x[6]^2 + 3*x[6]*x[7],
x[3]*x[7] - x[5]^2 + 2*x[5]*x[6] - 5*x[6]*x[7] + 3*x[7]^2,
x[4]*x[6] - x[5]^2 + 2*x[5]*x[7] + 2*x[6]^2 - 9*x[6]*x[7] + 6*x[7]^2,
x[4]*x[7] - x[5]*x[6] + 2*x[6]^2 - 4*x[6]*x[7] + 2*x[7]^2

w101 on X101 is given by:
Mapping from: Crv: X101 to Crv: X101
with equations : 
u[1]
u[2]
u[3]
u[4]
u[5]
u[6]
u[7]
-u[8]
and inverse
u[1]
u[2]
u[3]
u[4]
u[5]
u[6]
u[7]
-u[8]
Genus of X0(101) is  8
Genus of X0(101)/w101 is  1

We have found these points on X0(101):
[ (-1 : 0 : 0 : 0 : 0 : 0 : 0 : 1), (1 : 0 : 0 : 0 : 0 : 0 : 0 : 1) ]

Is Dtor := pts[1] - pts[2] a generator for J0(101)(Q)_tors?  true

X0(101)/w101 is actually the following elliptic curve:
Elliptic Curve defined by y^2 + 4*x*y + 7*y = x^3 - 10*x - 12 over Rational 
Field

It has the following MW group:
Abelian Group isomorphic to Z
Defined on 1 generator (free)

We have found  3  points on X_0(101)^2(Q).
1 of them are pullbacks of rationals from X0(101)/w101.
3
The number of possible cosets unknown points can land in is 2
5
The number of possible cosets unknown points can land in is 0
MWSieve achieved its goal? (true if succeeded, number otherwise)
true
It follows that if there is an unknown Q in X0(101)^2(Q), then 2[Q - bp] is 
fixed by w101.
This implies that [Q - bp] is fixed by w101 (since there is no 2-torsion)
Then Q ~ w101(Q), which implies that Q = w101(Q) because X0(101) is not 
hyperelliptic.
Then either Q is a pullback, or it is fixed by w101 pointwise.
If P = (X1:X2:X3:X4:X5:X6:X7:X8) is fixed by w101, then either X8 = 0 or P = 
(0:0:0:0:0:0:0:1) (not on X0(101))

We find all such P by imposing condition X8 = 0 and finding points on the 
scheme:

Scheme over Rational Field defined by
x[8],
x[1]^2 - 4*x[2]^2 - 4*x[2]*x[3] + 2*x[3]^2 + 8*x[3]*x[4] + 4*x[4]^2 - 
    20*x[4]*x[5] - 41*x[5]^2 + 112*x[5]*x[6] + 18*x[5]*x[7] - 22*x[6]^2 - 
    204*x[6]*x[7] + 170*x[7]^2 - x[8]^2,
x[1]*x[3] - x[2]^2 + 8*x[4]*x[5] - 8*x[5]^2 - 15*x[5]*x[6] + 21*x[5]*x[7] + 
    32*x[6]^2 - 61*x[6]*x[7] + 31*x[7]^2,
x[1]*x[4] - x[2]*x[3] + 8*x[4]*x[5] - 10*x[5]^2 - 5*x[5]*x[6] + 17*x[5]*x[7] + 
    22*x[6]^2 - 54*x[6]*x[7] + 30*x[7]^2,
x[1]*x[5] - x[3]^2 + 6*x[4]*x[5] - x[5]^2 - 20*x[5]*x[6] + 15*x[5]*x[7] + 
    25*x[6]^2 - 30*x[6]*x[7] + 9*x[7]^2,
x[1]*x[6] - x[3]*x[4] + 4*x[4]*x[5] - x[5]^2 - 10*x[5]*x[6] + 6*x[5]*x[7] + 
    12*x[6]^2 - 8*x[6]*x[7],
x[1]*x[7] - x[4]^2 + 2*x[4]*x[5] - x[5]*x[6] - x[5]*x[7] - 2*x[6]^2 + 
    10*x[6]*x[7] - 5*x[7]^2,
x[2]*x[4] - x[3]^2 + 8*x[5]^2 - 18*x[5]*x[6] - 2*x[5]*x[7] + 5*x[6]^2 + 
    26*x[6]*x[7] - 23*x[7]^2,
x[2]*x[5] - x[3]*x[4] + 6*x[5]^2 - 9*x[5]*x[6] - 6*x[5]*x[7] - 4*x[6]^2 + 
    33*x[6]*x[7] - 24*x[7]^2,
x[2]*x[6] - x[4]^2 + 4*x[5]^2 - x[5]*x[6] - 8*x[5]*x[7] - 10*x[6]^2 + 
    32*x[6]*x[7] - 20*x[7]^2,
x[2]*x[7] - x[4]*x[5] + 2*x[5]^2 - 3*x[5]*x[7] - 4*x[6]^2 + 12*x[6]*x[7] - 
    8*x[7]^2,
x[3]*x[5] - x[4]^2 + 6*x[5]*x[6] - 7*x[5]*x[7] - 10*x[6]^2 + 18*x[6]*x[7] - 
    8*x[7]^2,
x[3]*x[6] - x[4]*x[5] + 4*x[5]*x[6] - 2*x[5]*x[7] - 5*x[6]^2 + 3*x[6]*x[7],
x[3]*x[7] - x[5]^2 + 2*x[5]*x[6] - 5*x[6]*x[7] + 3*x[7]^2,
x[4]*x[6] - x[5]^2 + 2*x[5]*x[7] + 2*x[6]^2 - 9*x[6]*x[7] + 6*x[7]^2,
x[4]*x[7] - x[5]*x[6] + 2*x[6]^2 - 4*x[6]*x[7] + 2*x[7]^2

All such P are:
{@ (-312255/14117927*r1^13 + 7532589/14117927*r1^12 - 316552945/56471708*r1^11 +
    1928880373/56471708*r1^10 - 1919561629/14117927*r1^9 + 
    21292168973/56471708*r1^8 - 42582745737/56471708*r1^7 + 
    62429386129/56471708*r1^6 - 67479514287/56471708*r1^5 + 
    54248432759/56471708*r1^4 - 16901578683/28235854*r1^3 + 
    17753258123/56471708*r1^2 - 7785795003/56471708*r1 + 1705769353/56471708 : 
    136570/14117927*r1^13 - 3406016/14117927*r1^12 + 74874217/28235854*r1^11 - 
    241378129/14117927*r1^10 + 4104051993/56471708*r1^9 - 
    6118323011/28235854*r1^8 + 13242641357/28235854*r1^7 - 
    42298958117/56471708*r1^6 + 25003515291/28235854*r1^5 - 
    10886778599/14117927*r1^4 + 27935823331/56471708*r1^3 - 
    6774385543/28235854*r1^2 + 2554029659/28235854*r1 - 1123353655/56471708 : 
    110222/14117927*r1^13 - 2692677/14117927*r1^12 + 57923129/28235854*r1^11 - 
    733283731/56471708*r1^10 + 1544786437/28235854*r1^9 - 
    9281037935/56471708*r1^8 + 20669223075/56471708*r1^7 - 
    8692014556/14117927*r1^6 + 11069101751/14117927*r1^5 - 
    42034020153/56471708*r1^4 + 14269164887/28235854*r1^3 - 
    12522438083/56471708*r1^2 + 2734140353/56471708*r1 + 33044551/14117927 : 
    288448/14117927*r1^13 - 7180368/14117927*r1^12 + 78972663/14117927*r1^11 - 
    512583810/14117927*r1^10 + 8862323541/56471708*r1^9 - 
    6805761853/14117927*r1^8 + 61565480199/56471708*r1^7 - 
    104170691455/56471708*r1^6 + 33004611519/14117927*r1^5 - 
    61800927261/28235854*r1^4 + 82846115835/56471708*r1^3 - 
    9382403652/14117927*r1^2 + 10428667079/56471708*r1 - 1282484337/56471708 : 
    183919/14117927*r1^13 - 4493446/14117927*r1^12 + 193268761/56471708*r1^11 - 
    305418047/14117927*r1^10 + 1281178197/14117927*r1^9 - 
    15231801135/56471708*r1^8 + 8310356838/14117927*r1^7 - 
    54154147573/56471708*r1^6 + 65957694087/56471708*r1^5 - 
    14816535509/14117927*r1^4 + 19033312429/28235854*r1^3 - 
    16384070073/56471708*r1^2 + 1027691914/14117927*r1 - 271038965/56471708 : r1
: 1 : 0), (-312255/14117927*r2^13 + 7532589/14117927*r2^12 - 
    316552945/56471708*r2^11 + 1928880373/56471708*r2^10 - 
    1919561629/14117927*r2^9 + 21292168973/56471708*r2^8 - 
    42582745737/56471708*r2^7 + 62429386129/56471708*r2^6 - 
    67479514287/56471708*r2^5 + 54248432759/56471708*r2^4 - 
    16901578683/28235854*r2^3 + 17753258123/56471708*r2^2 - 
    7785795003/56471708*r2 + 1705769353/56471708 : 136570/14117927*r2^13 - 
    3406016/14117927*r2^12 + 74874217/28235854*r2^11 - 241378129/14117927*r2^10 
    + 4104051993/56471708*r2^9 - 6118323011/28235854*r2^8 + 
    13242641357/28235854*r2^7 - 42298958117/56471708*r2^6 + 
    25003515291/28235854*r2^5 - 10886778599/14117927*r2^4 + 
    27935823331/56471708*r2^3 - 6774385543/28235854*r2^2 + 
    2554029659/28235854*r2 - 1123353655/56471708 : 110222/14117927*r2^13 - 
    2692677/14117927*r2^12 + 57923129/28235854*r2^11 - 733283731/56471708*r2^10 
    + 1544786437/28235854*r2^9 - 9281037935/56471708*r2^8 + 
    20669223075/56471708*r2^7 - 8692014556/14117927*r2^6 + 
    11069101751/14117927*r2^5 - 42034020153/56471708*r2^4 + 
    14269164887/28235854*r2^3 - 12522438083/56471708*r2^2 + 
    2734140353/56471708*r2 + 33044551/14117927 : 288448/14117927*r2^13 - 
    7180368/14117927*r2^12 + 78972663/14117927*r2^11 - 512583810/14117927*r2^10 
    + 8862323541/56471708*r2^9 - 6805761853/14117927*r2^8 + 
    61565480199/56471708*r2^7 - 104170691455/56471708*r2^6 + 
    33004611519/14117927*r2^5 - 61800927261/28235854*r2^4 + 
    82846115835/56471708*r2^3 - 9382403652/14117927*r2^2 + 
    10428667079/56471708*r2 - 1282484337/56471708 : 183919/14117927*r2^13 - 
    4493446/14117927*r2^12 + 193268761/56471708*r2^11 - 305418047/14117927*r2^10
    + 1281178197/14117927*r2^9 - 15231801135/56471708*r2^8 + 
    8310356838/14117927*r2^7 - 54154147573/56471708*r2^6 + 
    65957694087/56471708*r2^5 - 14816535509/14117927*r2^4 + 
    19033312429/28235854*r2^3 - 16384070073/56471708*r2^2 + 
    1027691914/14117927*r2 - 271038965/56471708 : r2 : 1 : 0), 
(-312255/14117927*r3^13 + 7532589/14117927*r3^12 - 316552945/56471708*r3^11 + 
    1928880373/56471708*r3^10 - 1919561629/14117927*r3^9 + 
    21292168973/56471708*r3^8 - 42582745737/56471708*r3^7 + 
    62429386129/56471708*r3^6 - 67479514287/56471708*r3^5 + 
    54248432759/56471708*r3^4 - 16901578683/28235854*r3^3 + 
    17753258123/56471708*r3^2 - 7785795003/56471708*r3 + 1705769353/56471708 : 
    136570/14117927*r3^13 - 3406016/14117927*r3^12 + 74874217/28235854*r3^11 - 
    241378129/14117927*r3^10 + 4104051993/56471708*r3^9 - 
    6118323011/28235854*r3^8 + 13242641357/28235854*r3^7 - 
    42298958117/56471708*r3^6 + 25003515291/28235854*r3^5 - 
    10886778599/14117927*r3^4 + 27935823331/56471708*r3^3 - 
    6774385543/28235854*r3^2 + 2554029659/28235854*r3 - 1123353655/56471708 : 
    110222/14117927*r3^13 - 2692677/14117927*r3^12 + 57923129/28235854*r3^11 - 
    733283731/56471708*r3^10 + 1544786437/28235854*r3^9 - 
    9281037935/56471708*r3^8 + 20669223075/56471708*r3^7 - 
    8692014556/14117927*r3^6 + 11069101751/14117927*r3^5 - 
    42034020153/56471708*r3^4 + 14269164887/28235854*r3^3 - 
    12522438083/56471708*r3^2 + 2734140353/56471708*r3 + 33044551/14117927 : 
    288448/14117927*r3^13 - 7180368/14117927*r3^12 + 78972663/14117927*r3^11 - 
    512583810/14117927*r3^10 + 8862323541/56471708*r3^9 - 
    6805761853/14117927*r3^8 + 61565480199/56471708*r3^7 - 
    104170691455/56471708*r3^6 + 33004611519/14117927*r3^5 - 
    61800927261/28235854*r3^4 + 82846115835/56471708*r3^3 - 
    9382403652/14117927*r3^2 + 10428667079/56471708*r3 - 1282484337/56471708 : 
    183919/14117927*r3^13 - 4493446/14117927*r3^12 + 193268761/56471708*r3^11 - 
    305418047/14117927*r3^10 + 1281178197/14117927*r3^9 - 
    15231801135/56471708*r3^8 + 8310356838/14117927*r3^7 - 
    54154147573/56471708*r3^6 + 65957694087/56471708*r3^5 - 
    14816535509/14117927*r3^4 + 19033312429/28235854*r3^3 - 
    16384070073/56471708*r3^2 + 1027691914/14117927*r3 - 271038965/56471708 : r3
: 1 : 0), (-312255/14117927*r4^13 + 7532589/14117927*r4^12 - 
    316552945/56471708*r4^11 + 1928880373/56471708*r4^10 - 
    1919561629/14117927*r4^9 + 21292168973/56471708*r4^8 - 
    42582745737/56471708*r4^7 + 62429386129/56471708*r4^6 - 
    67479514287/56471708*r4^5 + 54248432759/56471708*r4^4 - 
    16901578683/28235854*r4^3 + 17753258123/56471708*r4^2 - 
    7785795003/56471708*r4 + 1705769353/56471708 : 136570/14117927*r4^13 - 
    3406016/14117927*r4^12 + 74874217/28235854*r4^11 - 241378129/14117927*r4^10 
    + 4104051993/56471708*r4^9 - 6118323011/28235854*r4^8 + 
    13242641357/28235854*r4^7 - 42298958117/56471708*r4^6 + 
    25003515291/28235854*r4^5 - 10886778599/14117927*r4^4 + 
    27935823331/56471708*r4^3 - 6774385543/28235854*r4^2 + 
    2554029659/28235854*r4 - 1123353655/56471708 : 110222/14117927*r4^13 - 
    2692677/14117927*r4^12 + 57923129/28235854*r4^11 - 733283731/56471708*r4^10 
    + 1544786437/28235854*r4^9 - 9281037935/56471708*r4^8 + 
    20669223075/56471708*r4^7 - 8692014556/14117927*r4^6 + 
    11069101751/14117927*r4^5 - 42034020153/56471708*r4^4 + 
    14269164887/28235854*r4^3 - 12522438083/56471708*r4^2 + 
    2734140353/56471708*r4 + 33044551/14117927 : 288448/14117927*r4^13 - 
    7180368/14117927*r4^12 + 78972663/14117927*r4^11 - 512583810/14117927*r4^10 
    + 8862323541/56471708*r4^9 - 6805761853/14117927*r4^8 + 
    61565480199/56471708*r4^7 - 104170691455/56471708*r4^6 + 
    33004611519/14117927*r4^5 - 61800927261/28235854*r4^4 + 
    82846115835/56471708*r4^3 - 9382403652/14117927*r4^2 + 
    10428667079/56471708*r4 - 1282484337/56471708 : 183919/14117927*r4^13 - 
    4493446/14117927*r4^12 + 193268761/56471708*r4^11 - 305418047/14117927*r4^10
    + 1281178197/14117927*r4^9 - 15231801135/56471708*r4^8 + 
    8310356838/14117927*r4^7 - 54154147573/56471708*r4^6 + 
    65957694087/56471708*r4^5 - 14816535509/14117927*r4^4 + 
    19033312429/28235854*r4^3 - 16384070073/56471708*r4^2 + 
    1027691914/14117927*r4 - 271038965/56471708 : r4 : 1 : 0), 
(-312255/14117927*r5^13 + 7532589/14117927*r5^12 - 316552945/56471708*r5^11 + 
    1928880373/56471708*r5^10 - 1919561629/14117927*r5^9 + 
    21292168973/56471708*r5^8 - 42582745737/56471708*r5^7 + 
    62429386129/56471708*r5^6 - 67479514287/56471708*r5^5 + 
    54248432759/56471708*r5^4 - 16901578683/28235854*r5^3 + 
    17753258123/56471708*r5^2 - 7785795003/56471708*r5 + 1705769353/56471708 : 
    136570/14117927*r5^13 - 3406016/14117927*r5^12 + 74874217/28235854*r5^11 - 
    241378129/14117927*r5^10 + 4104051993/56471708*r5^9 - 
    6118323011/28235854*r5^8 + 13242641357/28235854*r5^7 - 
    42298958117/56471708*r5^6 + 25003515291/28235854*r5^5 - 
    10886778599/14117927*r5^4 + 27935823331/56471708*r5^3 - 
    6774385543/28235854*r5^2 + 2554029659/28235854*r5 - 1123353655/56471708 : 
    110222/14117927*r5^13 - 2692677/14117927*r5^12 + 57923129/28235854*r5^11 - 
    733283731/56471708*r5^10 + 1544786437/28235854*r5^9 - 
    9281037935/56471708*r5^8 + 20669223075/56471708*r5^7 - 
    8692014556/14117927*r5^6 + 11069101751/14117927*r5^5 - 
    42034020153/56471708*r5^4 + 14269164887/28235854*r5^3 - 
    12522438083/56471708*r5^2 + 2734140353/56471708*r5 + 33044551/14117927 : 
    288448/14117927*r5^13 - 7180368/14117927*r5^12 + 78972663/14117927*r5^11 - 
    512583810/14117927*r5^10 + 8862323541/56471708*r5^9 - 
    6805761853/14117927*r5^8 + 61565480199/56471708*r5^7 - 
    104170691455/56471708*r5^6 + 33004611519/14117927*r5^5 - 
    61800927261/28235854*r5^4 + 82846115835/56471708*r5^3 - 
    9382403652/14117927*r5^2 + 10428667079/56471708*r5 - 1282484337/56471708 : 
    183919/14117927*r5^13 - 4493446/14117927*r5^12 + 193268761/56471708*r5^11 - 
    305418047/14117927*r5^10 + 1281178197/14117927*r5^9 - 
    15231801135/56471708*r5^8 + 8310356838/14117927*r5^7 - 
    54154147573/56471708*r5^6 + 65957694087/56471708*r5^5 - 
    14816535509/14117927*r5^4 + 19033312429/28235854*r5^3 - 
    16384070073/56471708*r5^2 + 1027691914/14117927*r5 - 271038965/56471708 : r5
: 1 : 0), (-312255/14117927*r6^13 + 7532589/14117927*r6^12 - 
    316552945/56471708*r6^11 + 1928880373/56471708*r6^10 - 
    1919561629/14117927*r6^9 + 21292168973/56471708*r6^8 - 
    42582745737/56471708*r6^7 + 62429386129/56471708*r6^6 - 
    67479514287/56471708*r6^5 + 54248432759/56471708*r6^4 - 
    16901578683/28235854*r6^3 + 17753258123/56471708*r6^2 - 
    7785795003/56471708*r6 + 1705769353/56471708 : 136570/14117927*r6^13 - 
    3406016/14117927*r6^12 + 74874217/28235854*r6^11 - 241378129/14117927*r6^10 
    + 4104051993/56471708*r6^9 - 6118323011/28235854*r6^8 + 
    13242641357/28235854*r6^7 - 42298958117/56471708*r6^6 + 
    25003515291/28235854*r6^5 - 10886778599/14117927*r6^4 + 
    27935823331/56471708*r6^3 - 6774385543/28235854*r6^2 + 
    2554029659/28235854*r6 - 1123353655/56471708 : 110222/14117927*r6^13 - 
    2692677/14117927*r6^12 + 57923129/28235854*r6^11 - 733283731/56471708*r6^10 
    + 1544786437/28235854*r6^9 - 9281037935/56471708*r6^8 + 
    20669223075/56471708*r6^7 - 8692014556/14117927*r6^6 + 
    11069101751/14117927*r6^5 - 42034020153/56471708*r6^4 + 
    14269164887/28235854*r6^3 - 12522438083/56471708*r6^2 + 
    2734140353/56471708*r6 + 33044551/14117927 : 288448/14117927*r6^13 - 
    7180368/14117927*r6^12 + 78972663/14117927*r6^11 - 512583810/14117927*r6^10 
    + 8862323541/56471708*r6^9 - 6805761853/14117927*r6^8 + 
    61565480199/56471708*r6^7 - 104170691455/56471708*r6^6 + 
    33004611519/14117927*r6^5 - 61800927261/28235854*r6^4 + 
    82846115835/56471708*r6^3 - 9382403652/14117927*r6^2 + 
    10428667079/56471708*r6 - 1282484337/56471708 : 183919/14117927*r6^13 - 
    4493446/14117927*r6^12 + 193268761/56471708*r6^11 - 305418047/14117927*r6^10
    + 1281178197/14117927*r6^9 - 15231801135/56471708*r6^8 + 
    8310356838/14117927*r6^7 - 54154147573/56471708*r6^6 + 
    65957694087/56471708*r6^5 - 14816535509/14117927*r6^4 + 
    19033312429/28235854*r6^3 - 16384070073/56471708*r6^2 + 
    1027691914/14117927*r6 - 271038965/56471708 : r6 : 1 : 0), 
(-312255/14117927*r7^13 + 7532589/14117927*r7^12 - 316552945/56471708*r7^11 + 
    1928880373/56471708*r7^10 - 1919561629/14117927*r7^9 + 
    21292168973/56471708*r7^8 - 42582745737/56471708*r7^7 + 
    62429386129/56471708*r7^6 - 67479514287/56471708*r7^5 + 
    54248432759/56471708*r7^4 - 16901578683/28235854*r7^3 + 
    17753258123/56471708*r7^2 - 7785795003/56471708*r7 + 1705769353/56471708 : 
    136570/14117927*r7^13 - 3406016/14117927*r7^12 + 74874217/28235854*r7^11 - 
    241378129/14117927*r7^10 + 4104051993/56471708*r7^9 - 
    6118323011/28235854*r7^8 + 13242641357/28235854*r7^7 - 
    42298958117/56471708*r7^6 + 25003515291/28235854*r7^5 - 
    10886778599/14117927*r7^4 + 27935823331/56471708*r7^3 - 
    6774385543/28235854*r7^2 + 2554029659/28235854*r7 - 1123353655/56471708 : 
    110222/14117927*r7^13 - 2692677/14117927*r7^12 + 57923129/28235854*r7^11 - 
    733283731/56471708*r7^10 + 1544786437/28235854*r7^9 - 
    9281037935/56471708*r7^8 + 20669223075/56471708*r7^7 - 
    8692014556/14117927*r7^6 + 11069101751/14117927*r7^5 - 
    42034020153/56471708*r7^4 + 14269164887/28235854*r7^3 - 
    12522438083/56471708*r7^2 + 2734140353/56471708*r7 + 33044551/14117927 : 
    288448/14117927*r7^13 - 7180368/14117927*r7^12 + 78972663/14117927*r7^11 - 
    512583810/14117927*r7^10 + 8862323541/56471708*r7^9 - 
    6805761853/14117927*r7^8 + 61565480199/56471708*r7^7 - 
    104170691455/56471708*r7^6 + 33004611519/14117927*r7^5 - 
    61800927261/28235854*r7^4 + 82846115835/56471708*r7^3 - 
    9382403652/14117927*r7^2 + 10428667079/56471708*r7 - 1282484337/56471708 : 
    183919/14117927*r7^13 - 4493446/14117927*r7^12 + 193268761/56471708*r7^11 - 
    305418047/14117927*r7^10 + 1281178197/14117927*r7^9 - 
    15231801135/56471708*r7^8 + 8310356838/14117927*r7^7 - 
    54154147573/56471708*r7^6 + 65957694087/56471708*r7^5 - 
    14816535509/14117927*r7^4 + 19033312429/28235854*r7^3 - 
    16384070073/56471708*r7^2 + 1027691914/14117927*r7 - 271038965/56471708 : r7
: 1 : 0), (-312255/14117927*r8^13 + 7532589/14117927*r8^12 - 
    316552945/56471708*r8^11 + 1928880373/56471708*r8^10 - 
    1919561629/14117927*r8^9 + 21292168973/56471708*r8^8 - 
    42582745737/56471708*r8^7 + 62429386129/56471708*r8^6 - 
    67479514287/56471708*r8^5 + 54248432759/56471708*r8^4 - 
    16901578683/28235854*r8^3 + 17753258123/56471708*r8^2 - 
    7785795003/56471708*r8 + 1705769353/56471708 : 136570/14117927*r8^13 - 
    3406016/14117927*r8^12 + 74874217/28235854*r8^11 - 241378129/14117927*r8^10 
    + 4104051993/56471708*r8^9 - 6118323011/28235854*r8^8 + 
    13242641357/28235854*r8^7 - 42298958117/56471708*r8^6 + 
    25003515291/28235854*r8^5 - 10886778599/14117927*r8^4 + 
    27935823331/56471708*r8^3 - 6774385543/28235854*r8^2 + 
    2554029659/28235854*r8 - 1123353655/56471708 : 110222/14117927*r8^13 - 
    2692677/14117927*r8^12 + 57923129/28235854*r8^11 - 733283731/56471708*r8^10 
    + 1544786437/28235854*r8^9 - 9281037935/56471708*r8^8 + 
    20669223075/56471708*r8^7 - 8692014556/14117927*r8^6 + 
    11069101751/14117927*r8^5 - 42034020153/56471708*r8^4 + 
    14269164887/28235854*r8^3 - 12522438083/56471708*r8^2 + 
    2734140353/56471708*r8 + 33044551/14117927 : 288448/14117927*r8^13 - 
    7180368/14117927*r8^12 + 78972663/14117927*r8^11 - 512583810/14117927*r8^10 
    + 8862323541/56471708*r8^9 - 6805761853/14117927*r8^8 + 
    61565480199/56471708*r8^7 - 104170691455/56471708*r8^6 + 
    33004611519/14117927*r8^5 - 61800927261/28235854*r8^4 + 
    82846115835/56471708*r8^3 - 9382403652/14117927*r8^2 + 
    10428667079/56471708*r8 - 1282484337/56471708 : 183919/14117927*r8^13 - 
    4493446/14117927*r8^12 + 193268761/56471708*r8^11 - 305418047/14117927*r8^10
    + 1281178197/14117927*r8^9 - 15231801135/56471708*r8^8 + 
    8310356838/14117927*r8^7 - 54154147573/56471708*r8^6 + 
    65957694087/56471708*r8^5 - 14816535509/14117927*r8^4 + 
    19033312429/28235854*r8^3 - 16384070073/56471708*r8^2 + 
    1027691914/14117927*r8 - 271038965/56471708 : r8 : 1 : 0), 
(-312255/14117927*r9^13 + 7532589/14117927*r9^12 - 316552945/56471708*r9^11 + 
    1928880373/56471708*r9^10 - 1919561629/14117927*r9^9 + 
    21292168973/56471708*r9^8 - 42582745737/56471708*r9^7 + 
    62429386129/56471708*r9^6 - 67479514287/56471708*r9^5 + 
    54248432759/56471708*r9^4 - 16901578683/28235854*r9^3 + 
    17753258123/56471708*r9^2 - 7785795003/56471708*r9 + 1705769353/56471708 : 
    136570/14117927*r9^13 - 3406016/14117927*r9^12 + 74874217/28235854*r9^11 - 
    241378129/14117927*r9^10 + 4104051993/56471708*r9^9 - 
    6118323011/28235854*r9^8 + 13242641357/28235854*r9^7 - 
    42298958117/56471708*r9^6 + 25003515291/28235854*r9^5 - 
    10886778599/14117927*r9^4 + 27935823331/56471708*r9^3 - 
    6774385543/28235854*r9^2 + 2554029659/28235854*r9 - 1123353655/56471708 : 
    110222/14117927*r9^13 - 2692677/14117927*r9^12 + 57923129/28235854*r9^11 - 
    733283731/56471708*r9^10 + 1544786437/28235854*r9^9 - 
    9281037935/56471708*r9^8 + 20669223075/56471708*r9^7 - 
    8692014556/14117927*r9^6 + 11069101751/14117927*r9^5 - 
    42034020153/56471708*r9^4 + 14269164887/28235854*r9^3 - 
    12522438083/56471708*r9^2 + 2734140353/56471708*r9 + 33044551/14117927 : 
    288448/14117927*r9^13 - 7180368/14117927*r9^12 + 78972663/14117927*r9^11 - 
    512583810/14117927*r9^10 + 8862323541/56471708*r9^9 - 
    6805761853/14117927*r9^8 + 61565480199/56471708*r9^7 - 
    104170691455/56471708*r9^6 + 33004611519/14117927*r9^5 - 
    61800927261/28235854*r9^4 + 82846115835/56471708*r9^3 - 
    9382403652/14117927*r9^2 + 10428667079/56471708*r9 - 1282484337/56471708 : 
    183919/14117927*r9^13 - 4493446/14117927*r9^12 + 193268761/56471708*r9^11 - 
    305418047/14117927*r9^10 + 1281178197/14117927*r9^9 - 
    15231801135/56471708*r9^8 + 8310356838/14117927*r9^7 - 
    54154147573/56471708*r9^6 + 65957694087/56471708*r9^5 - 
    14816535509/14117927*r9^4 + 19033312429/28235854*r9^3 - 
    16384070073/56471708*r9^2 + 1027691914/14117927*r9 - 271038965/56471708 : r9
: 1 : 0), (-312255/14117927*r10^13 + 7532589/14117927*r10^12 - 
    316552945/56471708*r10^11 + 1928880373/56471708*r10^10 - 
    1919561629/14117927*r10^9 + 21292168973/56471708*r10^8 - 
    42582745737/56471708*r10^7 + 62429386129/56471708*r10^6 - 
    67479514287/56471708*r10^5 + 54248432759/56471708*r10^4 - 
    16901578683/28235854*r10^3 + 17753258123/56471708*r10^2 - 
    7785795003/56471708*r10 + 1705769353/56471708 : 136570/14117927*r10^13 - 
    3406016/14117927*r10^12 + 74874217/28235854*r10^11 - 
    241378129/14117927*r10^10 + 4104051993/56471708*r10^9 - 
    6118323011/28235854*r10^8 + 13242641357/28235854*r10^7 - 
    42298958117/56471708*r10^6 + 25003515291/28235854*r10^5 - 
    10886778599/14117927*r10^4 + 27935823331/56471708*r10^3 - 
    6774385543/28235854*r10^2 + 2554029659/28235854*r10 - 1123353655/56471708 : 
110222/14117927*r10^13 - 2692677/14117927*r10^12 + 57923129/28235854*r10^11 - 
    733283731/56471708*r10^10 + 1544786437/28235854*r10^9 - 
    9281037935/56471708*r10^8 + 20669223075/56471708*r10^7 - 
    8692014556/14117927*r10^6 + 11069101751/14117927*r10^5 - 
    42034020153/56471708*r10^4 + 14269164887/28235854*r10^3 - 
    12522438083/56471708*r10^2 + 2734140353/56471708*r10 + 33044551/14117927 : 
    288448/14117927*r10^13 - 7180368/14117927*r10^12 + 78972663/14117927*r10^11 
    - 512583810/14117927*r10^10 + 8862323541/56471708*r10^9 - 
    6805761853/14117927*r10^8 + 61565480199/56471708*r10^7 - 
    104170691455/56471708*r10^6 + 33004611519/14117927*r10^5 - 
    61800927261/28235854*r10^4 + 82846115835/56471708*r10^3 - 
    9382403652/14117927*r10^2 + 10428667079/56471708*r10 - 1282484337/56471708 :
183919/14117927*r10^13 - 4493446/14117927*r10^12 + 193268761/56471708*r10^11 - 
    305418047/14117927*r10^10 + 1281178197/14117927*r10^9 - 
    15231801135/56471708*r10^8 + 8310356838/14117927*r10^7 - 
    54154147573/56471708*r10^6 + 65957694087/56471708*r10^5 - 
    14816535509/14117927*r10^4 + 19033312429/28235854*r10^3 - 
    16384070073/56471708*r10^2 + 1027691914/14117927*r10 - 271038965/56471708 : 
r10 : 1 : 0), (-312255/14117927*r11^13 + 7532589/14117927*r11^12 - 
    316552945/56471708*r11^11 + 1928880373/56471708*r11^10 - 
    1919561629/14117927*r11^9 + 21292168973/56471708*r11^8 - 
    42582745737/56471708*r11^7 + 62429386129/56471708*r11^6 - 
    67479514287/56471708*r11^5 + 54248432759/56471708*r11^4 - 
    16901578683/28235854*r11^3 + 17753258123/56471708*r11^2 - 
    7785795003/56471708*r11 + 1705769353/56471708 : 136570/14117927*r11^13 - 
    3406016/14117927*r11^12 + 74874217/28235854*r11^11 - 
    241378129/14117927*r11^10 + 4104051993/56471708*r11^9 - 
    6118323011/28235854*r11^8 + 13242641357/28235854*r11^7 - 
    42298958117/56471708*r11^6 + 25003515291/28235854*r11^5 - 
    10886778599/14117927*r11^4 + 27935823331/56471708*r11^3 - 
    6774385543/28235854*r11^2 + 2554029659/28235854*r11 - 1123353655/56471708 : 
110222/14117927*r11^13 - 2692677/14117927*r11^12 + 57923129/28235854*r11^11 - 
    733283731/56471708*r11^10 + 1544786437/28235854*r11^9 - 
    9281037935/56471708*r11^8 + 20669223075/56471708*r11^7 - 
    8692014556/14117927*r11^6 + 11069101751/14117927*r11^5 - 
    42034020153/56471708*r11^4 + 14269164887/28235854*r11^3 - 
    12522438083/56471708*r11^2 + 2734140353/56471708*r11 + 33044551/14117927 : 
    288448/14117927*r11^13 - 7180368/14117927*r11^12 + 78972663/14117927*r11^11 
    - 512583810/14117927*r11^10 + 8862323541/56471708*r11^9 - 
    6805761853/14117927*r11^8 + 61565480199/56471708*r11^7 - 
    104170691455/56471708*r11^6 + 33004611519/14117927*r11^5 - 
    61800927261/28235854*r11^4 + 82846115835/56471708*r11^3 - 
    9382403652/14117927*r11^2 + 10428667079/56471708*r11 - 1282484337/56471708 :
183919/14117927*r11^13 - 4493446/14117927*r11^12 + 193268761/56471708*r11^11 - 
    305418047/14117927*r11^10 + 1281178197/14117927*r11^9 - 
    15231801135/56471708*r11^8 + 8310356838/14117927*r11^7 - 
    54154147573/56471708*r11^6 + 65957694087/56471708*r11^5 - 
    14816535509/14117927*r11^4 + 19033312429/28235854*r11^3 - 
    16384070073/56471708*r11^2 + 1027691914/14117927*r11 - 271038965/56471708 : 
r11 : 1 : 0), (-312255/14117927*r12^13 + 7532589/14117927*r12^12 - 
    316552945/56471708*r12^11 + 1928880373/56471708*r12^10 - 
    1919561629/14117927*r12^9 + 21292168973/56471708*r12^8 - 
    42582745737/56471708*r12^7 + 62429386129/56471708*r12^6 - 
    67479514287/56471708*r12^5 + 54248432759/56471708*r12^4 - 
    16901578683/28235854*r12^3 + 17753258123/56471708*r12^2 - 
    7785795003/56471708*r12 + 1705769353/56471708 : 136570/14117927*r12^13 - 
    3406016/14117927*r12^12 + 74874217/28235854*r12^11 - 
    241378129/14117927*r12^10 + 4104051993/56471708*r12^9 - 
    6118323011/28235854*r12^8 + 13242641357/28235854*r12^7 - 
    42298958117/56471708*r12^6 + 25003515291/28235854*r12^5 - 
    10886778599/14117927*r12^4 + 27935823331/56471708*r12^3 - 
    6774385543/28235854*r12^2 + 2554029659/28235854*r12 - 1123353655/56471708 : 
110222/14117927*r12^13 - 2692677/14117927*r12^12 + 57923129/28235854*r12^11 - 
    733283731/56471708*r12^10 + 1544786437/28235854*r12^9 - 
    9281037935/56471708*r12^8 + 20669223075/56471708*r12^7 - 
    8692014556/14117927*r12^6 + 11069101751/14117927*r12^5 - 
    42034020153/56471708*r12^4 + 14269164887/28235854*r12^3 - 
    12522438083/56471708*r12^2 + 2734140353/56471708*r12 + 33044551/14117927 : 
    288448/14117927*r12^13 - 7180368/14117927*r12^12 + 78972663/14117927*r12^11 
    - 512583810/14117927*r12^10 + 8862323541/56471708*r12^9 - 
    6805761853/14117927*r12^8 + 61565480199/56471708*r12^7 - 
    104170691455/56471708*r12^6 + 33004611519/14117927*r12^5 - 
    61800927261/28235854*r12^4 + 82846115835/56471708*r12^3 - 
    9382403652/14117927*r12^2 + 10428667079/56471708*r12 - 1282484337/56471708 :
183919/14117927*r12^13 - 4493446/14117927*r12^12 + 193268761/56471708*r12^11 - 
    305418047/14117927*r12^10 + 1281178197/14117927*r12^9 - 
    15231801135/56471708*r12^8 + 8310356838/14117927*r12^7 - 
    54154147573/56471708*r12^6 + 65957694087/56471708*r12^5 - 
    14816535509/14117927*r12^4 + 19033312429/28235854*r12^3 - 
    16384070073/56471708*r12^2 + 1027691914/14117927*r12 - 271038965/56471708 : 
r12 : 1 : 0), (-312255/14117927*r13^13 + 7532589/14117927*r13^12 - 
    316552945/56471708*r13^11 + 1928880373/56471708*r13^10 - 
    1919561629/14117927*r13^9 + 21292168973/56471708*r13^8 - 
    42582745737/56471708*r13^7 + 62429386129/56471708*r13^6 - 
    67479514287/56471708*r13^5 + 54248432759/56471708*r13^4 - 
    16901578683/28235854*r13^3 + 17753258123/56471708*r13^2 - 
    7785795003/56471708*r13 + 1705769353/56471708 : 136570/14117927*r13^13 - 
    3406016/14117927*r13^12 + 74874217/28235854*r13^11 - 
    241378129/14117927*r13^10 + 4104051993/56471708*r13^9 - 
    6118323011/28235854*r13^8 + 13242641357/28235854*r13^7 - 
    42298958117/56471708*r13^6 + 25003515291/28235854*r13^5 - 
    10886778599/14117927*r13^4 + 27935823331/56471708*r13^3 - 
    6774385543/28235854*r13^2 + 2554029659/28235854*r13 - 1123353655/56471708 : 
110222/14117927*r13^13 - 2692677/14117927*r13^12 + 57923129/28235854*r13^11 - 
    733283731/56471708*r13^10 + 1544786437/28235854*r13^9 - 
    9281037935/56471708*r13^8 + 20669223075/56471708*r13^7 - 
    8692014556/14117927*r13^6 + 11069101751/14117927*r13^5 - 
    42034020153/56471708*r13^4 + 14269164887/28235854*r13^3 - 
    12522438083/56471708*r13^2 + 2734140353/56471708*r13 + 33044551/14117927 : 
    288448/14117927*r13^13 - 7180368/14117927*r13^12 + 78972663/14117927*r13^11 
    - 512583810/14117927*r13^10 + 8862323541/56471708*r13^9 - 
    6805761853/14117927*r13^8 + 61565480199/56471708*r13^7 - 
    104170691455/56471708*r13^6 + 33004611519/14117927*r13^5 - 
    61800927261/28235854*r13^4 + 82846115835/56471708*r13^3 - 
    9382403652/14117927*r13^2 + 10428667079/56471708*r13 - 1282484337/56471708 :
183919/14117927*r13^13 - 4493446/14117927*r13^12 + 193268761/56471708*r13^11 - 
    305418047/14117927*r13^10 + 1281178197/14117927*r13^9 - 
    15231801135/56471708*r13^8 + 8310356838/14117927*r13^7 - 
    54154147573/56471708*r13^6 + 65957694087/56471708*r13^5 - 
    14816535509/14117927*r13^4 + 19033312429/28235854*r13^3 - 
    16384070073/56471708*r13^2 + 1027691914/14117927*r13 - 271038965/56471708 : 
r13 : 1 : 0), (-312255/14117927*r14^13 + 7532589/14117927*r14^12 - 
    316552945/56471708*r14^11 + 1928880373/56471708*r14^10 - 
    1919561629/14117927*r14^9 + 21292168973/56471708*r14^8 - 
    42582745737/56471708*r14^7 + 62429386129/56471708*r14^6 - 
    67479514287/56471708*r14^5 + 54248432759/56471708*r14^4 - 
    16901578683/28235854*r14^3 + 17753258123/56471708*r14^2 - 
    7785795003/56471708*r14 + 1705769353/56471708 : 136570/14117927*r14^13 - 
    3406016/14117927*r14^12 + 74874217/28235854*r14^11 - 
    241378129/14117927*r14^10 + 4104051993/56471708*r14^9 - 
    6118323011/28235854*r14^8 + 13242641357/28235854*r14^7 - 
    42298958117/56471708*r14^6 + 25003515291/28235854*r14^5 - 
    10886778599/14117927*r14^4 + 27935823331/56471708*r14^3 - 
    6774385543/28235854*r14^2 + 2554029659/28235854*r14 - 1123353655/56471708 : 
110222/14117927*r14^13 - 2692677/14117927*r14^12 + 57923129/28235854*r14^11 - 
    733283731/56471708*r14^10 + 1544786437/28235854*r14^9 - 
    9281037935/56471708*r14^8 + 20669223075/56471708*r14^7 - 
    8692014556/14117927*r14^6 + 11069101751/14117927*r14^5 - 
    42034020153/56471708*r14^4 + 14269164887/28235854*r14^3 - 
    12522438083/56471708*r14^2 + 2734140353/56471708*r14 + 33044551/14117927 : 
    288448/14117927*r14^13 - 7180368/14117927*r14^12 + 78972663/14117927*r14^11 
    - 512583810/14117927*r14^10 + 8862323541/56471708*r14^9 - 
    6805761853/14117927*r14^8 + 61565480199/56471708*r14^7 - 
    104170691455/56471708*r14^6 + 33004611519/14117927*r14^5 - 
    61800927261/28235854*r14^4 + 82846115835/56471708*r14^3 - 
    9382403652/14117927*r14^2 + 10428667079/56471708*r14 - 1282484337/56471708 :
183919/14117927*r14^13 - 4493446/14117927*r14^12 + 193268761/56471708*r14^11 - 
    305418047/14117927*r14^10 + 1281178197/14117927*r14^9 - 
    15231801135/56471708*r14^8 + 8310356838/14117927*r14^7 - 
    54154147573/56471708*r14^6 + 65957694087/56471708*r14^5 - 
    14816535509/14117927*r14^4 + 19033312429/28235854*r14^3 - 
    16384070073/56471708*r14^2 + 1027691914/14117927*r14 - 271038965/56471708 : 
r14 : 1 : 0) @}

Clearly, none such points are quadratic.
Hence there are no quadratic points on X0(101) not coming from pullbacks of 
rationals.

Total time: 246.419 seconds, Total memory usage: 160.22MB
