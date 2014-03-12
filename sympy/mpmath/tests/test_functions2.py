import math
from sympy.mpmath import *

def test_bessel():
    mp.dps = 15
    assert j0(1).ae(0.765197686557966551)
    assert j0(pi).ae(-0.304242177644093864)
    assert j0(1000).ae(0.0247866861524201746)
    assert j0(-25).ae(0.0962667832759581162)
    assert j1(1).ae(0.440050585744933516)
    assert j1(pi).ae(0.284615343179752757)
    assert j1(1000).ae(0.00472831190708952392)
    assert j1(-25).ae(0.125350249580289905)
    assert besselj(5,1).ae(0.000249757730211234431)
    assert besselj(5+0j,1).ae(0.000249757730211234431)
    assert besselj(5,pi).ae(0.0521411843671184747)
    assert besselj(5,1000).ae(0.00502540694523318607)
    assert besselj(5,-25).ae(0.0660079953984229934)
    assert besselj(-3,2).ae(-0.128943249474402051)
    assert besselj(-4,2).ae(0.0339957198075684341)
    assert besselj(3,3+2j).ae(0.424718794929639595942 + 0.625665327745785804812j)
    assert besselj(0.25,4).ae(-0.374760630804249715)
    assert besselj(1+2j,3+4j).ae(0.319247428741872131 - 0.669557748880365678j)
    assert (besselj(3, 10**10) * 10**5).ae(0.76765081748139204023)
    assert bessely(-0.5, 0) == 0
    assert bessely(0.5, 0) == -inf
    assert bessely(1.5, 0) == -inf
    assert bessely(0,0) == -inf
    assert bessely(-0.4, 0) == -inf
    assert bessely(-0.6, 0) == inf
    assert bessely(-1, 0) == inf
    assert bessely(-1.4, 0) == inf
    assert bessely(-1.6, 0) == -inf
    assert bessely(-1, 0) == inf
    assert bessely(-2, 0) == -inf
    assert bessely(-3, 0) == inf
    assert bessely(0.5, 0) == -inf
    assert bessely(1, 0) == -inf
    assert bessely(1.5, 0) == -inf
    assert bessely(2, 0) == -inf
    assert bessely(2.5, 0) == -inf
    assert bessely(3, 0) == -inf
    assert bessely(0,0.5).ae(-0.44451873350670655715)
    assert bessely(1,0.5).ae(-1.4714723926702430692)
    assert bessely(-1,0.5).ae(1.4714723926702430692)
    assert bessely(3.5,0.5).ae(-138.86400867242488443)
    assert bessely(0,3+4j).ae(4.6047596915010138655-8.8110771408232264208j)
    assert bessely(0,j).ae(-0.26803248203398854876+1.26606587775200833560j)
    assert (bessely(3, 10**10) * 10**5).ae(0.21755917537013204058)
    assert besseli(0,0) == 1
    assert besseli(1,0) == 0
    assert besseli(2,0) == 0
    assert besseli(-1,0) == 0
    assert besseli(-2,0) == 0
    assert besseli(0,0.5).ae(1.0634833707413235193)
    assert besseli(1,0.5).ae(0.25789430539089631636)
    assert besseli(-1,0.5).ae(0.25789430539089631636)
    assert besseli(3.5,0.5).ae(0.00068103597085793815863)
    assert besseli(0,3+4j).ae(-3.3924877882755196097-1.3239458916287264815j)
    assert besseli(0,j).ae(besselj(0,1))
    assert (besseli(3, 10**10) * mpf(10)**(-4342944813)).ae(4.2996028505491271875)
    assert besselk(0,0) == inf
    assert besselk(1,0) == inf
    assert besselk(2,0) == inf
    assert besselk(-1,0) == inf
    assert besselk(-2,0) == inf
    assert besselk(0,0.5).ae(0.92441907122766586178)
    assert besselk(1,0.5).ae(1.6564411200033008937)
    assert besselk(-1,0.5).ae(1.6564411200033008937)
    assert besselk(3.5,0.5).ae(207.48418747548460607)
    assert besselk(0,3+4j).ae(-0.007239051213570155013+0.026510418350267677215j)
    assert besselk(0,j).ae(-0.13863371520405399968-1.20196971531720649914j)
    assert (besselk(3, 10**10) * mpf(10)**4342944824).ae(1.1628981033356187851)

def test_bessel_zeros():
    mp.dps = 15
    assert besseljzero(0,1).ae(2.40482555769577276869)
    assert besseljzero(2,1).ae(5.1356223018406825563)
    assert besseljzero(1,50).ae(157.86265540193029781)
    assert besseljzero(10,1).ae(14.475500686554541220)
    assert besseljzero(0.5,3).ae(9.4247779607693797153)
    assert besseljzero(2,1,1).ae(3.0542369282271403228)
    assert besselyzero(0,1).ae(0.89357696627916752158)
    assert besselyzero(2,1).ae(3.3842417671495934727)
    assert besselyzero(1,50).ae(156.29183520147840108)
    assert besselyzero(10,1).ae(12.128927704415439387)
    assert besselyzero(0.5,3).ae(7.8539816339744830962)
    assert besselyzero(2,1,1).ae(5.0025829314460639452)

def test_hankel():
    mp.dps = 15
    assert hankel1(0,0.5).ae(0.93846980724081290423-0.44451873350670655715j)
    assert hankel1(1,0.5).ae(0.2422684576748738864-1.4714723926702430692j)
    assert hankel1(-1,0.5).ae(-0.2422684576748738864+1.4714723926702430692j)
    assert hankel1(1.5,0.5).ae(0.0917016996256513026-2.5214655504213378514j)
    assert hankel1(1.5,3+4j).ae(0.0066806866476728165382-0.0036684231610839127106j)
    assert hankel2(0,0.5).ae(0.93846980724081290423+0.44451873350670655715j)
    assert hankel2(1,0.5).ae(0.2422684576748738864+1.4714723926702430692j)
    assert hankel2(-1,0.5).ae(-0.2422684576748738864-1.4714723926702430692j)
    assert hankel2(1.5,0.5).ae(0.0917016996256513026+2.5214655504213378514j)
    assert hankel2(1.5,3+4j).ae(14.783528526098567526-7.397390270853446512j)

def test_struve():
    mp.dps = 15
    assert struveh(2,3).ae(0.74238666967748318564)
    assert struveh(-2.5,3).ae(0.41271003220971599344)
    assert struvel(2,3).ae(1.7476573277362782744)
    assert struvel(-2.5,3).ae(1.5153394466819651377)

def test_whittaker():
    mp.dps = 15
    assert whitm(2,3,4).ae(49.753745589025246591)
    assert whitw(2,3,4).ae(14.111656223052932215)

def test_kelvin():
    mp.dps = 15
    assert ber(2,3).ae(0.80836846563726819091)
    assert ber(3,4).ae(-0.28262680167242600233)
    assert ber(-3,2).ae(-0.085611448496796363669)
    assert bei(2,3).ae(-0.89102236377977331571)
    assert bei(-3,2).ae(-0.14420994155731828415)
    assert ker(2,3).ae(0.12839126695733458928)
    assert ker(-3,2).ae(-0.29802153400559142783)
    assert ker(0.5,3).ae(-0.085662378535217097524)
    assert kei(2,3).ae(0.036804426134164634000)
    assert kei(-3,2).ae(0.88682069845786731114)
    assert kei(0.5,3).ae(0.013633041571314302948)

def test_hyper_misc():
    mp.dps = 15
    assert hyp0f1(1,0) == 1
    assert hyp1f1(1,2,0) == 1
    assert hyp1f2(1,2,3,0) == 1
    assert hyp2f1(1,2,3,0) == 1
    assert hyp2f2(1,2,3,4,0) == 1
    assert hyp2f3(1,2,3,4,5,0) == 1
    # Degenerate case: 0F0
    assert hyper([],[],0) == 1
    assert hyper([],[],-2).ae(exp(-2))
    # Degenerate case: 1F0
    assert hyper([2],[],1.5) == 4
    #
    assert hyp2f1((1,3),(2,3),(5,6),mpf(27)/32).ae(1.6)
    assert hyp2f1((1,4),(1,2),(3,4),mpf(80)/81).ae(1.8)
    assert hyp2f1((2,3),(1,1),(3,2),(2+j)/3).ae(1.327531603558679093+0.439585080092769253j)
    mp.dps = 25
    v = mpc('1.2282306665029814734863026', '-0.1225033830118305184672133')
    assert hyper([(3,4),2+j,1],[1,5,j/3],mpf(1)/5+j/8).ae(v)
    mp.dps = 15

def test_elliptic_integrals():
    mp.dps = 15
    assert ellipk(0).ae(pi/2)
    assert ellipk(0.5).ae(gamma(0.25)**2/(4*sqrt(pi)))
    assert ellipk(1) == inf
    assert ellipk(1+0j) == inf
    assert ellipk(-1).ae('1.3110287771460599052')
    assert ellipk(-2).ae('1.1714200841467698589')
    assert isinstance(ellipk(-2), mpf)
    assert isinstance(ellipe(-2), mpf)
    assert ellipk(-50).ae('0.47103424540873331679')
    mp.dps = 30
    n1 = +fraction(99999,100000)
    n2 = +fraction(100001,100000)
    mp.dps = 15
    assert ellipk(n1).ae('7.1427724505817781901')
    assert ellipk(n2).ae(mpc('7.1427417367963090109', '-1.5707923998261688019'))
    assert ellipe(n1).ae('1.0000332138990829170')
    v = ellipe(n2)
    assert v.real.ae('0.999966786328145474069137')
    assert (v.imag*10**6).ae('7.853952181727432')
    assert ellipk(2).ae(mpc('1.3110287771460599052', '-1.3110287771460599052'))
    assert ellipk(50).ae(mpc('0.22326753950210985451', '-0.47434723226254522087'))
    assert ellipk(3+4j).ae(mpc('0.91119556380496500866', '0.63133428324134524388'))
    assert ellipk(3-4j).ae(mpc('0.91119556380496500866', '-0.63133428324134524388'))
    assert ellipk(-3+4j).ae(mpc('0.95357894880405122483', '0.23093044503746114444'))
    assert ellipk(-3-4j).ae(mpc('0.95357894880405122483', '-0.23093044503746114444'))
    assert isnan(ellipk(nan))
    assert isnan(ellipe(nan))
    assert ellipk(inf) == 0
    assert isinstance(ellipk(inf), mpc)
    assert ellipk(-inf) == 0
    assert ellipk(1+0j) == inf
    assert ellipe(0).ae(pi/2)
    assert ellipe(0.5).ae(pi**(mpf(3)/2)/gamma(0.25)**2 +gamma(0.25)**2/(8*sqrt(pi)))
    assert ellipe(1) == 1
    assert ellipe(1+0j) == 1
    assert ellipe(inf) == mpc(0,inf)
    assert ellipe(-inf) == inf
    assert ellipe(3+4j).ae(1.4995535209333469543-1.5778790079127582745j)
    assert ellipe(3-4j).ae(1.4995535209333469543+1.5778790079127582745j)
    assert ellipe(-3+4j).ae(2.5804237855343377803-0.8306096791000413778j)
    assert ellipe(-3-4j).ae(2.5804237855343377803+0.8306096791000413778j)
    assert ellipe(2).ae(0.59907011736779610372+0.59907011736779610372j)
    assert ellipe('1e-1000000000').ae(pi/2)
    assert ellipk('1e-1000000000').ae(pi/2)
    assert ellipe(-pi).ae(2.4535865983838923)
    mp.dps = 50
    assert ellipk(1/pi).ae('1.724756270009501831744438120951614673874904182624739673')
    assert ellipe(1/pi).ae('1.437129808135123030101542922290970050337425479058225712')
    assert ellipk(-10*pi).ae('0.5519067523886233967683646782286965823151896970015484512')
    assert ellipe(-10*pi).ae('5.926192483740483797854383268707108012328213431657645509')
    v = ellipk(pi)
    assert v.real.ae('0.973089521698042334840454592642137667227167622330325225')
    assert v.imag.ae('-1.156151296372835303836814390793087600271609993858798016')
    v = ellipe(pi)
    assert v.real.ae('0.4632848917264710404078033487934663562998345622611263332')
    assert v.imag.ae('1.0637961621753130852473300451583414489944099504180510966')
    mp.dps = 15

def test_exp_integrals():
    mp.dps = 15
    x = +e
    z = e + sqrt(3)*j
    assert ei(x).ae(8.21168165538361560)
    assert li(x).ae(1.89511781635593676)
    assert si(x).ae(1.82104026914756705)
    assert ci(x).ae(0.213958001340379779)
    assert shi(x).ae(4.11520706247846193)
    assert chi(x).ae(4.09647459290515367)
    assert fresnels(x).ae(0.437189718149787643)
    assert fresnelc(x).ae(0.401777759590243012)
    assert airyai(x).ae(0.0108502401568586681)
    assert airybi(x).ae(8.98245748585468627)
    assert ei(z).ae(3.72597969491314951 + 7.34213212314224421j)
    assert li(z).ae(2.28662658112562502 + 1.50427225297269364j)
    assert si(z).ae(2.48122029237669054 + 0.12684703275254834j)
    assert ci(z).ae(0.169255590269456633 - 0.892020751420780353j)
    assert shi(z).ae(1.85810366559344468 + 3.66435842914920263j)
    assert chi(z).ae(1.86787602931970484 + 3.67777369399304159j)
    assert fresnels(z/3).ae(0.034534397197008182 + 0.754859844188218737j)
    assert fresnelc(z/3).ae(1.261581645990027372 + 0.417949198775061893j)
    assert airyai(z).ae(-0.0162552579839056062 - 0.0018045715700210556j)
    assert airybi(z).ae(-4.98856113282883371 + 2.08558537872180623j)
    assert li(0) == 0.0
    assert li(1) == -inf
    assert li(inf) == inf
    assert isinstance(li(0.7), mpf)
    assert si(inf).ae(pi/2)
    assert si(-inf).ae(-pi/2)
    assert ci(inf) == 0
    assert ci(0) == -inf
    assert isinstance(ei(-0.7), mpf)
    assert airyai(inf) == 0
    assert airybi(inf) == inf
    assert airyai(-inf) == 0
    assert airybi(-inf) == 0
    assert fresnels(inf) == 0.5
    assert fresnelc(inf) == 0.5
    assert fresnels(-inf) == -0.5
    assert fresnelc(-inf) == -0.5
    assert shi(0) == 0
    assert shi(inf) == inf
    assert shi(-inf) == -inf
    assert chi(0) == -inf
    assert chi(inf) == inf

def test_ei():
    mp.dps = 15
    assert ei(0) == -inf
    assert ei(inf) == inf
    assert ei(-inf) == -0.0
    assert ei(20+70j).ae(6.1041351911152984397e6 - 2.7324109310519928872e6j)
    # tests for the asymptotic expansion
    # values checked with Mathematica ExpIntegralEi
    mp.dps = 50
    r = ei(20000)
    s = '3.8781962825045010930273870085501819470698476975019e+8681'
    assert str(r) == s
    r = ei(-200)
    s = '-6.8852261063076355977108174824557929738368086933303e-90'
    assert str(r) == s
    r =ei(20000 + 10*j)
    sre = '-3.255138234032069402493850638874410725961401274106e+8681'
    sim = '-2.1081929993474403520785942429469187647767369645423e+8681'
    assert str(r.real) == sre and str(r.imag) == sim
    mp.dps = 15
    # More asymptotic expansions
    assert chi(-10**6+100j).ae('1.3077239389562548386e+434288 + 7.6808956999707408158e+434287j')
    assert shi(-10**6+100j).ae('-1.3077239389562548386e+434288 - 7.6808956999707408158e+434287j')
    mp.dps = 15
    assert ei(10j).ae(-0.0454564330044553726+3.2291439210137706686j)
    assert ei(100j).ae(-0.0051488251426104921+3.1330217936839529126j)
    u = ei(fmul(10**20, j, exact=True))
    assert u.real.ae(-6.4525128526578084421345e-21, abs_eps=0, rel_eps=8*eps)
    assert u.imag.ae(pi)
    assert ei(-10j).ae(-0.0454564330044553726-3.2291439210137706686j)
    assert ei(-100j).ae(-0.0051488251426104921-3.1330217936839529126j)
    u = ei(fmul(-10**20, j, exact=True))
    assert u.real.ae(-6.4525128526578084421345e-21, abs_eps=0, rel_eps=8*eps)
    assert u.imag.ae(-pi)
    assert ei(10+10j).ae(-1576.1504265768517448+436.9192317011328140j)
    u = ei(-10+10j)
    assert u.real.ae(7.6698978415553488362543e-7, abs_eps=0, rel_eps=8*eps)
    assert u.imag.ae(3.141595611735621062025)

def test_e1():
    mp.dps = 15
    assert e1(0) == inf
    assert e1(inf) == 0
    assert e1(-inf) == mpc(-inf, -pi)
    assert e1(10j).ae(0.045456433004455372635 + 0.087551267423977430100j)
    assert e1(100j).ae(0.0051488251426104921444 - 0.0085708599058403258790j)
    assert e1(fmul(10**20, j, exact=True)).ae(6.4525128526578084421e-21 - 7.6397040444172830039e-21j, abs_eps=0, rel_eps=8*eps)
    assert e1(-10j).ae(0.045456433004455372635 - 0.087551267423977430100j)
    assert e1(-100j).ae(0.0051488251426104921444 + 0.0085708599058403258790j)
    assert e1(fmul(-10**20, j, exact=True)).ae(6.4525128526578084421e-21 + 7.6397040444172830039e-21j, abs_eps=0, rel_eps=8*eps)

def test_expint():
    mp.dps = 15
    assert expint(0,0) == inf
    assert expint(0,1).ae(1/e)
    assert expint(0,1.5).ae(2/exp(1.5)/3)
    assert expint(1,1).ae(-ei(-1))
    assert expint(2,0).ae(1)
    assert expint(3,0).ae(1/2.)
    assert expint(4,0).ae(1/3.)
    assert expint(-2, 0.5).ae(26/sqrt(e))
    assert expint(-1,-1) == 0
    assert expint(-2,-1).ae(-e)
    assert expint(5.5, 0).ae(2/9.)
    assert expint(2.00000001,0).ae(100000000./100000001)
    assert expint(2+3j,4-j).ae(0.0023461179581675065414+0.0020395540604713669262j)
    assert expint('1.01', '1e-1000').ae(99.9999999899412802)
    assert expint('1.000000000001', 3.5).ae(0.00697013985754701819446)
    assert expint(2,3).ae(3*ei(-3)+exp(-3))
    assert (expint(10,20)*10**10).ae(0.694439055541231353)
    assert expint(3,inf) == 0
    assert expint(3.2,inf) == 0
    assert expint(3.2+2j,inf) == 0
    assert expint(1,3j).ae(-0.11962978600800032763 + 0.27785620120457163717j)
    assert expint(1,3).ae(0.013048381094197037413)
    assert expint(1,-3).ae(-ei(3)-pi*j)
    #assert expint(3) == expint(1,3)
    assert expint(1,-20).ae(-25615652.66405658882 - 3.1415926535897932385j)
    assert expint(1000000,0).ae(1./999999)
    assert expint(0,2+3j).ae(-0.025019798357114678171 + 0.027980439405104419040j)
    assert expint(-1,2+3j).ae(-0.022411973626262070419 + 0.038058922011377716932j)
    assert expint(-1.5,0) == inf

def test_trig_integrals():
    mp.dps = 30
    assert si(mpf(1)/1000000).ae('0.000000999999999999944444444444446111')
    assert ci(mpf(1)/1000000).ae('-13.2382948930629912435014366276')
    assert si(10**10).ae('1.5707963267075846569685111517747537')
    assert ci(10**10).ae('-4.87506025174822653785729773959e-11')
    assert si(10**100).ae(pi/2)
    assert (ci(10**100)*10**100).ae('-0.372376123661276688262086695553')
    assert si(-3) == -si(3)
    assert ci(-3).ae(ci(3) + pi*j)
    # Test complex structure
    mp.dps = 15
    assert mp.ci(50).ae(-0.0056283863241163054402)
    assert mp.ci(50+2j).ae(-0.018378282946133067149+0.070352808023688336193j)
    assert mp.ci(20j).ae(1.28078263320282943611e7+1.5707963267949j)
    assert mp.ci(-2+20j).ae(-4.050116856873293505e6+1.207476188206989909e7j)
    assert mp.ci(-50+2j).ae(-0.0183782829461330671+3.0712398455661049023j)
    assert mp.ci(-50).ae(-0.0056283863241163054+3.1415926535897932385j)
    assert mp.ci(-50-2j).ae(-0.0183782829461330671-3.0712398455661049023j)
    assert mp.ci(-2-20j).ae(-4.050116856873293505e6-1.207476188206989909e7j)
    assert mp.ci(-20j).ae(1.28078263320282943611e7-1.5707963267949j)
    assert mp.ci(50-2j).ae(-0.018378282946133067149-0.070352808023688336193j)
    assert mp.si(50).ae(1.5516170724859358947)
    assert mp.si(50+2j).ae(1.497884414277228461-0.017515007378437448j)
    assert mp.si(20j).ae(1.2807826332028294459e7j)
    assert mp.si(-2+20j).ae(-1.20747603112735722103e7-4.050116856873293554e6j)
    assert mp.si(-50+2j).ae(-1.497884414277228461-0.017515007378437448j)
    assert mp.si(-50).ae(-1.5516170724859358947)
    assert mp.si(-50-2j).ae(-1.497884414277228461+0.017515007378437448j)
    assert mp.si(-2-20j).ae(-1.20747603112735722103e7+4.050116856873293554e6j)
    assert mp.si(-20j).ae(-1.2807826332028294459e7j)
    assert mp.si(50-2j).ae(1.497884414277228461+0.017515007378437448j)
    assert mp.chi(50j).ae(-0.0056283863241163054+1.5707963267948966192j)
    assert mp.chi(-2+50j).ae(-0.0183782829461330671+1.6411491348185849554j)
    assert mp.chi(-20).ae(1.28078263320282943611e7+3.1415926535898j)
    assert mp.chi(-20-2j).ae(-4.050116856873293505e6+1.20747571696809187053e7j)
    assert mp.chi(-2-50j).ae(-0.0183782829461330671-1.6411491348185849554j)
    assert mp.chi(-50j).ae(-0.0056283863241163054-1.5707963267948966192j)
    assert mp.chi(2-50j).ae(-0.0183782829461330671-1.500443518771208283j)
    assert mp.chi(20-2j).ae(-4.050116856873293505e6-1.20747603112735722951e7j)
    assert mp.chi(20).ae(1.2807826332028294361e7)
    assert mp.chi(2+50j).ae(-0.0183782829461330671+1.500443518771208283j)
    assert mp.shi(50j).ae(1.5516170724859358947j)
    assert mp.shi(-2+50j).ae(0.017515007378437448+1.497884414277228461j)
    assert mp.shi(-20).ae(-1.2807826332028294459e7)
    assert mp.shi(-20-2j).ae(4.050116856873293554e6-1.20747603112735722103e7j)
    assert mp.shi(-2-50j).ae(0.017515007378437448-1.497884414277228461j)
    assert mp.shi(-50j).ae(-1.5516170724859358947j)
    assert mp.shi(2-50j).ae(-0.017515007378437448-1.497884414277228461j)
    assert mp.shi(20-2j).ae(-4.050116856873293554e6-1.20747603112735722103e7j)
    assert mp.shi(20).ae(1.2807826332028294459e7)
    assert mp.shi(2+50j).ae(-0.017515007378437448+1.497884414277228461j)
    def ae(x,y,tol=1e-12):
        return abs(x-y) <= abs(y)*tol
    assert fp.ci(fp.inf) == 0
    assert ae(fp.ci(fp.ninf), fp.pi*1j)
    assert ae(fp.si(fp.inf), fp.pi/2)
    assert ae(fp.si(fp.ninf), -fp.pi/2)
    assert fp.si(0) == 0
    assert ae(fp.ci(50), -0.0056283863241163054402)
    assert ae(fp.ci(50+2j), -0.018378282946133067149+0.070352808023688336193j)
    assert ae(fp.ci(20j), 1.28078263320282943611e7+1.5707963267949j)
    assert ae(fp.ci(-2+20j), -4.050116856873293505e6+1.207476188206989909e7j)
    assert ae(fp.ci(-50+2j), -0.0183782829461330671+3.0712398455661049023j)
    assert ae(fp.ci(-50), -0.0056283863241163054+3.1415926535897932385j)
    assert ae(fp.ci(-50-2j), -0.0183782829461330671-3.0712398455661049023j)
    assert ae(fp.ci(-2-20j), -4.050116856873293505e6-1.207476188206989909e7j)
    assert ae(fp.ci(-20j), 1.28078263320282943611e7-1.5707963267949j)
    assert ae(fp.ci(50-2j), -0.018378282946133067149-0.070352808023688336193j)
    assert ae(fp.si(50), 1.5516170724859358947)
    assert ae(fp.si(50+2j), 1.497884414277228461-0.017515007378437448j)
    assert ae(fp.si(20j), 1.2807826332028294459e7j)
    assert ae(fp.si(-2+20j), -1.20747603112735722103e7-4.050116856873293554e6j)
    assert ae(fp.si(-50+2j), -1.497884414277228461-0.017515007378437448j)
    assert ae(fp.si(-50), -1.5516170724859358947)
    assert ae(fp.si(-50-2j), -1.497884414277228461+0.017515007378437448j)
    assert ae(fp.si(-2-20j), -1.20747603112735722103e7+4.050116856873293554e6j)
    assert ae(fp.si(-20j), -1.2807826332028294459e7j)
    assert ae(fp.si(50-2j), 1.497884414277228461+0.017515007378437448j)
    assert ae(fp.chi(50j), -0.0056283863241163054+1.5707963267948966192j)
    assert ae(fp.chi(-2+50j), -0.0183782829461330671+1.6411491348185849554j)
    assert ae(fp.chi(-20), 1.28078263320282943611e7+3.1415926535898j)
    assert ae(fp.chi(-20-2j), -4.050116856873293505e6+1.20747571696809187053e7j)
    assert ae(fp.chi(-2-50j), -0.0183782829461330671-1.6411491348185849554j)
    assert ae(fp.chi(-50j), -0.0056283863241163054-1.5707963267948966192j)
    assert ae(fp.chi(2-50j), -0.0183782829461330671-1.500443518771208283j)
    assert ae(fp.chi(20-2j), -4.050116856873293505e6-1.20747603112735722951e7j)
    assert ae(fp.chi(20), 1.2807826332028294361e7)
    assert ae(fp.chi(2+50j), -0.0183782829461330671+1.500443518771208283j)
    assert ae(fp.shi(50j), 1.5516170724859358947j)
    assert ae(fp.shi(-2+50j), 0.017515007378437448+1.497884414277228461j)
    assert ae(fp.shi(-20), -1.2807826332028294459e7)
    assert ae(fp.shi(-20-2j), 4.050116856873293554e6-1.20747603112735722103e7j)
    assert ae(fp.shi(-2-50j), 0.017515007378437448-1.497884414277228461j)
    assert ae(fp.shi(-50j), -1.5516170724859358947j)
    assert ae(fp.shi(2-50j), -0.017515007378437448-1.497884414277228461j)
    assert ae(fp.shi(20-2j), -4.050116856873293554e6-1.20747603112735722103e7j)
    assert ae(fp.shi(20), 1.2807826332028294459e7)
    assert ae(fp.shi(2+50j), -0.017515007378437448+1.497884414277228461j)

def test_airy():
    mp.dps = 15
    assert (airyai(10)*10**10).ae(1.1047532552898687)
    assert (airybi(10)/10**9).ae(0.45564115354822515)
    assert (airyai(1000)*10**9158).ae(9.306933063179556004)
    assert (airybi(1000)/10**9154).ae(5.4077118391949465477)
    assert airyai(-1000).ae(0.055971895773019918842)
    assert airybi(-1000).ae(-0.083264574117080633012)
    assert (airyai(100+100j)*10**188).ae(2.9099582462207032076 + 2.353013591706178756j)
    assert (airybi(100+100j)/10**185).ae(1.7086751714463652039 - 3.1416590020830804578j)

def test_hyper_0f1():
    mp.dps = 15
    v = 8.63911136507950465
    assert hyper([],[(1,3)],1.5).ae(v)
    assert hyper([],[1/3.],1.5).ae(v)
    assert hyp0f1(1/3.,1.5).ae(v)
    assert hyp0f1((1,3),1.5).ae(v)
    # Asymptotic expansion
    assert hyp0f1(3,1e9).ae('4.9679055380347771271e+27455')
    assert hyp0f1(3,1e9j).ae('-2.1222788784457702157e+19410 + 5.0840597555401854116e+19410j')

def test_hyper_1f1():
    mp.dps = 15
    v = 1.2917526488617656673
    assert hyper([(1,2)],[(3,2)],0.7).ae(v)
    assert hyper([(1,2)],[(3,2)],0.7+0j).ae(v)
    assert hyper([0.5],[(3,2)],0.7).ae(v)
    assert hyper([0.5],[1.5],0.7).ae(v)
    assert hyper([0.5],[(3,2)],0.7+0j).ae(v)
    assert hyper([0.5],[1.5],0.7+0j).ae(v)
    assert hyper([(1,2)],[1.5+0j],0.7).ae(v)
    assert hyper([0.5+0j],[1.5],0.7).ae(v)
    assert hyper([0.5+0j],[1.5+0j],0.7+0j).ae(v)
    assert hyp1f1(0.5,1.5,0.7).ae(v)
    assert hyp1f1((1,2),1.5,0.7).ae(v)
    # Asymptotic expansion
    assert hyp1f1(2,3,1e10).ae('2.1555012157015796988e+4342944809')
    assert (hyp1f1(2,3,1e10j)*10**10).ae(-0.97501205020039745852 - 1.7462392454512132074j)
    # Shouldn't use asymptotic expansion
    assert hyp1f1(-2, 1, 10000).ae(49980001)
    # Bug
    assert hyp1f1(1j,fraction(1,3),0.415-69.739j).ae(25.857588206024346592 + 15.738060264515292063j)

def test_hyper_2f1():
    mp.dps = 15
    v = 1.0652207633823291032
    assert hyper([(1,2), (3,4)], [2], 0.3).ae(v)
    assert hyper([(1,2), 0.75], [2], 0.3).ae(v)
    assert hyper([0.5, 0.75], [2.0], 0.3).ae(v)
    assert hyper([0.5, 0.75], [2.0], 0.3+0j).ae(v)
    assert hyper([0.5+0j, (3,4)], [2.0], 0.3+0j).ae(v)
    assert hyper([0.5+0j, (3,4)], [2.0], 0.3).ae(v)
    assert hyper([0.5, (3,4)], [2.0+0j], 0.3).ae(v)
    assert hyper([0.5+0j, 0.75+0j], [2.0+0j], 0.3+0j).ae(v)
    v = 1.09234681096223231717 + 0.18104859169479360380j
    assert hyper([(1,2),0.75+j], [2], 0.5).ae(v)
    assert hyper([0.5,0.75+j], [2.0], 0.5).ae(v)
    assert hyper([0.5,0.75+j], [2.0], 0.5+0j).ae(v)
    assert hyper([0.5,0.75+j], [2.0+0j], 0.5+0j).ae(v)
    v = 0.9625 - 0.125j
    assert hyper([(3,2),-1],[4], 0.1+j/3).ae(v)
    assert hyper([1.5,-1.0],[4], 0.1+j/3).ae(v)
    assert hyper([1.5,-1.0],[4+0j], 0.1+j/3).ae(v)
    assert hyper([1.5+0j,-1.0+0j],[4+0j], 0.1+j/3).ae(v)
    v = 1.02111069501693445001 - 0.50402252613466859521j
    assert hyper([(2,10),(3,10)],[(4,10)],1.5).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],1.5).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],1.5+0j).ae(v)
    v = 0.76922501362865848528 + 0.32640579593235886194j
    assert hyper([(2,10),(3,10)],[(4,10)],4+2j).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],4+2j).ae(v)
    assert hyper([0.2,(3,10)],[(4,10)],4+2j).ae(v)

def test_hyper_2f1_hard():
    mp.dps = 15
    # Singular cases
    assert hyp2f1(2,-1,-1,3).ae(0.25)
    assert hyp2f1(2,-2,-2,3).ae(0.25)
    assert hyp2f1(2,-1,-1,3,eliminate=False) == 7
    assert hyp2f1(2,-2,-2,3,eliminate=False) == 34
    assert hyp2f1(2,-2,-3,3) == 14
    assert hyp2f1(2,-3,-2,3) == inf
    assert hyp2f1(2,-1.5,-1.5,3) == 0.25
    assert hyp2f1(1,2,3,0) == 1
    assert hyp2f1(0,1,0,0) == 1
    assert hyp2f1(0,0,0,0) == 1
    assert isnan(hyp2f1(1,1,0,0))
    assert hyp2f1(2,-1,-5, 0.25+0.25j).ae(1.1+0.1j)
    assert hyp2f1(2,-5,-5, 0.25+0.25j, eliminate=False).ae(163./128 + 125./128*j)
    assert hyp2f1(0.7235, -1, -5, 0.3).ae(1.04341)
    assert hyp2f1(0.7235, -5, -5, 0.3, eliminate=False).ae(1.2939225017815903812)
    assert hyp2f1(-1,-2,4,1) == 1.5
    assert hyp2f1(1,2,-3,1) == inf
    assert hyp2f1(-2,-2,1,1) == 6
    assert hyp2f1(1,-2,-4,1).ae(5./3)
    assert hyp2f1(0,-6,-4,1) == 1
    assert hyp2f1(0,-3,-4,1) == 1
    assert hyp2f1(0,0,0,1) == 1
    assert hyp2f1(1,0,0,1,eliminate=False) == 1
    assert hyp2f1(1,1,0,1) == inf
    assert hyp2f1(1,-6,-4,1) == inf
    assert hyp2f1(-7.2,-0.5,-4.5,1) == 0
    assert hyp2f1(-7.2,-1,-2,1).ae(-2.6)
    assert hyp2f1(1,-0.5,-4.5, 1) == inf
    assert hyp2f1(1,0.5,-4.5, 1) == -inf
    # Check evaluation on / close to unit circle
    z = exp(j*pi/3)
    w = (nthroot(2,3)+1)*exp(j*pi/12)/nthroot(3,4)**3
    assert hyp2f1('1/2','1/6','1/3', z).ae(w)
    assert hyp2f1('1/2','1/6','1/3', z.conjugate()).ae(w.conjugate())
    assert hyp2f1(0.25, (1,3), 2, '0.999').ae(1.06826449496030635)
    assert hyp2f1(0.25, (1,3), 2, '1.001').ae(1.06867299254830309446-0.00001446586793975874j)
    assert hyp2f1(0.25, (1,3), 2, -1).ae(0.96656584492524351673)
    assert hyp2f1(0.25, (1,3), 2, j).ae(0.99041766248982072266+0.03777135604180735522j)
    assert hyp2f1(2,3,5,'0.99').ae(27.699347904322690602)
    assert hyp2f1((3,2),-0.5,3,'0.99').ae(0.68403036843911661388)
    assert hyp2f1(2,3,5,1j).ae(0.37290667145974386127+0.59210004902748285917j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), 0.95*exp(j*k)) for k in range(1,15)]).ae(52.851400204289452922+6.244285013912953225j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), 1.05*exp(j*k)) for k in range(1,15)]).ae(54.506013786220655330-3.000118813413217097j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), exp(j*k)) for k in range(1,15)]).ae(55.792077935955314887+1.731986485778500241j)
    assert hyp2f1(2,2.5,-3.25,0.999).ae(218373932801217082543180041.33)
    # Branches
    assert hyp2f1(1,1,2,1.01).ae(4.5595744415723676911-3.1104877758314784539j)
    assert hyp2f1(1,1,2,1.01+0.1j).ae(2.4149427480552782484+1.4148224796836938829j)
    assert hyp2f1(1,1,2,3+4j).ae(0.14576709331407297807+0.48379185417980360773j)
    assert hyp2f1(1,1,2,4).ae(-0.27465307216702742285 - 0.78539816339744830962j)
    assert hyp2f1(1,1,2,-4).ae(0.40235947810852509365)
    # Other:
    # Cancellation with a large parameter involved (bug reported on sage-devel)
    assert hyp2f1(112, (51,10), (-9,10), -0.99999).ae(-1.6241361047970862961e-24, abs_eps=0, rel_eps=eps*16)

def test_hyper_3f2_etc():
    assert hyper([1,2,3],[1.5,8],-1).ae(0.67108992351533333030)
    assert hyper([1,2,3,4],[5,6,7], -1).ae(0.90232988035425506008)
    assert hyper([1,2,3],[1.25,5], 1).ae(28.924181329701905701)
    assert hyper([1,2,3,4],[5,6,7],5).ae(1.5192307344006649499-1.1529845225075537461j)
    assert hyper([1,2,3,4,5],[6,7,8,9],-1).ae(0.96288759462882357253)
    assert hyper([1,2,3,4,5],[6,7,8,9],1).ae(1.0428697385885855841)
    assert hyper([1,2,3,4,5],[6,7,8,9],5).ae(1.33980653631074769423-0.07143405251029226699j)
    assert hyper([1,2.79,3.08,4.37],[5.2,6.1,7.3],5).ae(1.0996321464692607231-1.7748052293979985001j)
    assert hyper([1,1,1],[1,2],1) == inf
    assert hyper([1,1,1],[2,(101,100)],1).ae(100.01621213528313220)
    # slow -- covered by doctests
    #assert hyper([1,1,1],[2,3],0.9999).ae(1.2897972005319693905)

def test_hyper_u():
    mp.dps = 15
    assert hyperu(2,-3,0).ae(0.05)
    assert hyperu(2,-3.5,0).ae(4./99)
    assert hyperu(2,0,0) == 0.5
    assert hyperu(-5,1,0) == -120
    assert hyperu(-5,2,0) == inf
    assert hyperu(-5,-2,0) == 0
    assert hyperu(7,7,3).ae(0.00014681269365593503986)  #exp(3)*gammainc(-6,3)
    assert hyperu(2,-3,4).ae(0.011836478100271995559)
    assert hyperu(3,4,5).ae(1./125)
    assert hyperu(2,3,0.0625) == 256
    assert hyperu(-1,2,0.25+0.5j) == -1.75+0.5j
    assert hyperu(0.5,1.5,7.25).ae(2/sqrt(29))
    assert hyperu(2,6,pi).ae(0.55804439825913399130)
    assert (hyperu((3,2),8,100+201j)*10**4).ae(-0.3797318333856738798 - 2.9974928453561707782j)
    assert (hyperu((5,2),(-1,2),-5000)*10**10).ae(-5.6681877926881664678j)
    # XXX: fails because of undetected cancellation in low level series code
    # Alternatively: could use asymptotic series here, if convergence test
    # tweaked back to recognize this one
    #assert (hyperu((5,2),(-1,2),-500)*10**7).ae(-1.82526906001593252847j)

def test_hyper_2f0():
    mp.dps = 15
    assert hyper([1,2],[],3) == hyp2f0(1,2,3)
    assert hyp2f0(2,3,7).ae(0.0116108068639728714668 - 0.0073727413865865802130j)
    assert hyp2f0(2,3,0) == 1
    assert hyp2f0(0,0,0) == 1
    assert hyp2f0(-1,-1,1).ae(2)
    assert hyp2f0(-4,1,1.5).ae(62.5)
    assert hyp2f0(-4,1,50).ae(147029801)
    assert hyp2f0(-4,1,0.0001).ae(0.99960011997600240000)
    assert hyp2f0(0.5,0.25,0.001).ae(1.0001251174078538115)
    assert hyp2f0(0.5,0.25,3+4j).ae(0.85548875824755163518 + 0.21636041283392292973j)
    # Important: cancellation check
    assert hyp2f0((1,6),(5,6),-0.02371708245126284498).ae(0.996785723120804309)
    # Should be exact; polynomial case
    assert hyp2f0(-2,1,0.5+0.5j,zeroprec=200) == 0
    assert hyp2f0(1,-2,0.5+0.5j,zeroprec=200) == 0
    # There used to be a bug in thresholds that made one of the following hang
    for d in [15, 50, 80]:
        mp.dps = d
        assert hyp2f0(1.5, 0.5, 0.009).ae('1.006867007239309717945323585695344927904000945829843527398772456281301440034218290443367270629519483 + 1.238277162240704919639384945859073461954721356062919829456053965502443570466701567100438048602352623e-46j')

def test_hyper_1f2():
    mp.dps = 15
    assert hyper([1],[2,3],4) == hyp1f2(1,2,3,4)
    a1,b1,b2 = (1,10),(2,3),1./16
    assert hyp1f2(a1,b1,b2,10).ae(298.7482725554557568)
    assert hyp1f2(a1,b1,b2,100).ae(224128961.48602947604)
    assert hyp1f2(a1,b1,b2,1000).ae(1.1669528298622675109e+27)
    assert hyp1f2(a1,b1,b2,10000).ae(2.4780514622487212192e+86)
    assert hyp1f2(a1,b1,b2,100000).ae(1.3885391458871523997e+274)
    assert hyp1f2(a1,b1,b2,1000000).ae('9.8851796978960318255e+867')
    assert hyp1f2(a1,b1,b2,10**7).ae('1.1505659189516303646e+2746')
    assert hyp1f2(a1,b1,b2,10**8).ae('1.4672005404314334081e+8685')
    assert hyp1f2(a1,b1,b2,10**20).ae('3.6888217332150976493e+8685889636')
    assert hyp1f2(a1,b1,b2,10*j).ae(-16.163252524618572878 - 44.321567896480184312j)
    assert hyp1f2(a1,b1,b2,100*j).ae(61938.155294517848171 + 637349.45215942348739j)
    assert hyp1f2(a1,b1,b2,1000*j).ae(8455057657257695958.7 + 6261969266997571510.6j)
    assert hyp1f2(a1,b1,b2,10000*j).ae(-8.9771211184008593089e+60 + 4.6550528111731631456e+59j)
    assert hyp1f2(a1,b1,b2,100000*j).ae(2.6398091437239324225e+193 + 4.1658080666870618332e+193j)
    assert hyp1f2(a1,b1,b2,1000000*j).ae('3.5999042951925965458e+613 + 1.5026014707128947992e+613j')
    assert hyp1f2(a1,b1,b2,10**7*j).ae('-8.3208715051623234801e+1939 - 3.6752883490851869429e+1941j')
    assert hyp1f2(a1,b1,b2,10**8*j).ae('2.0724195707891484454e+6140 - 1.3276619482724266387e+6141j')
    assert hyp1f2(a1,b1,b2,10**20*j).ae('-1.1734497974795488504e+6141851462 + 1.1498106965385471542e+6141851462j')

def test_hyper_2f3():
    mp.dps = 15
    assert hyper([1,2],[3,4,5],6) == hyp2f3(1,2,3,4,5,6)
    a1,a2,b1,b2,b3 = (1,10),(2,3),(3,10), 2, 1./16
    # Check asymptotic expansion
    assert hyp2f3(a1,a2,b1,b2,b3,10).ae(128.98207160698659976)
    assert hyp2f3(a1,a2,b1,b2,b3,1000).ae(6.6309632883131273141e25)
    assert hyp2f3(a1,a2,b1,b2,b3,10000).ae(4.6863639362713340539e84)
    assert hyp2f3(a1,a2,b1,b2,b3,100000).ae(8.6632451236103084119e271)
    assert hyp2f3(a1,a2,b1,b2,b3,10**6).ae('2.0291718386574980641e865')
    assert hyp2f3(a1,a2,b1,b2,b3,10**7).ae('7.7639836665710030977e2742')
    assert hyp2f3(a1,a2,b1,b2,b3,10**8).ae('3.2537462584071268759e8681')
    assert hyp2f3(a1,a2,b1,b2,b3,10**20).ae('1.2966030542911614163e+8685889627')
    assert hyp2f3(a1,a2,b1,b2,b3,10*j).ae(-18.551602185587547854 - 13.348031097874113552j)
    assert hyp2f3(a1,a2,b1,b2,b3,100*j).ae(78634.359124504488695 + 74459.535945281973996j)
    assert hyp2f3(a1,a2,b1,b2,b3,1000*j).ae(597682550276527901.59 - 65136194809352613.078j)
    assert hyp2f3(a1,a2,b1,b2,b3,10000*j).ae(-1.1779696326238582496e+59 + 1.2297607505213133872e+59j)
    assert hyp2f3(a1,a2,b1,b2,b3,100000*j).ae(2.9844228969804380301e+191 + 7.5587163231490273296e+190j)
    assert hyp2f3(a1,a2,b1,b2,b3,1000000*j).ae('7.4859161049322370311e+610 - 2.8467477015940090189e+610j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**7*j).ae('-1.7477645579418800826e+1938 - 1.7606522995808116405e+1938j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**8*j).ae('-1.6932731942958401784e+6137 - 2.4521909113114629368e+6137j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**20*j).ae('-2.0988815677627225449e+6141851451 + 5.7708223542739208681e+6141851452j')

def test_hyper_2f2():
    mp.dps = 15
    assert hyper([1,2],[3,4],5) == hyp2f2(1,2,3,4,5)
    a1,a2,b1,b2 = (3,10),4,(1,2),1./16
    assert hyp2f2(a1,a2,b1,b2,10).ae(448225936.3377556696)
    assert hyp2f2(a1,a2,b1,b2,10000).ae('1.2012553712966636711e+4358')
    assert hyp2f2(a1,a2,b1,b2,-20000).ae(-0.04182343755661214626)
    assert hyp2f2(a1,a2,b1,b2,10**20).ae('1.1148680024303263661e+43429448190325182840')

def test_orthpoly():
    mp.dps = 15
    assert jacobi(-4,2,3,0.7).ae(22800./4913)
    assert jacobi(3,2,4,5.5) == 4133.125
    assert jacobi(1.5,5/6.,4,0).ae(-1.0851951434075508417)
    assert jacobi(-2, 1, 2, 4).ae(-0.16)
    assert jacobi(2, -1, 2.5, 4).ae(34.59375)
    #assert jacobi(2, -1, 2, 4) == 28.5
    assert legendre(5, 7) == 129367
    assert legendre(0.5,0).ae(0.53935260118837935667)
    assert legendre(-1,-1) == 1
    assert legendre(0,-1) == 1
    assert legendre(0, 1) == 1
    assert legendre(1, -1) == -1
    assert legendre(7, 1) == 1
    assert legendre(7, -1) == -1
    assert legendre(8,1.5).ae(15457523./32768)
    assert legendre(j,-j).ae(2.4448182735671431011 + 0.6928881737669934843j)
    assert chebyu(5,1) == 6
    assert chebyt(3,2) == 26
    assert legendre(3.5,-1) == inf
    assert legendre(4.5,-1) == -inf
    assert legendre(3.5+1j,-1) == mpc(inf,inf)
    assert legendre(4.5+1j,-1) == mpc(-inf,-inf)
    assert laguerre(4, -2, 3).ae(-1.125)
    assert laguerre(3, 1+j, 0.5).ae(0.2291666666666666667 + 2.5416666666666666667j)

def test_hermite():
    mp.dps = 15
    assert hermite(-2, 0).ae(0.5)
    assert hermite(-1, 0).ae(0.88622692545275801365)
    assert hermite(0, 0).ae(1)
    assert hermite(1, 0) == 0
    assert hermite(2, 0).ae(-2)
    assert hermite(0, 2).ae(1)
    assert hermite(1, 2).ae(4)
    assert hermite(1, -2).ae(-4)
    assert hermite(2, -2).ae(14)
    assert hermite(0.5, 0).ae(0.69136733903629335053)
    assert hermite(9, 0) == 0
    assert hermite(4,4).ae(3340)
    assert hermite(3,4).ae(464)
    assert hermite(-4,4).ae(0.00018623860287512396181)
    assert hermite(-3,4).ae(0.0016540169879668766270)
    assert hermite(9, 2.5j).ae(13638725j)
    assert hermite(9, -2.5j).ae(-13638725j)
    assert hermite(9, 100).ae(511078883759363024000)
    assert hermite(9, -100).ae(-511078883759363024000)
    assert hermite(9, 100j).ae(512922083920643024000j)
    assert hermite(9, -100j).ae(-512922083920643024000j)
    assert hermite(-9.5, 2.5j).ae(-2.9004951258126778174e-6 + 1.7601372934039951100e-6j)
    assert hermite(-9.5, -2.5j).ae(-2.9004951258126778174e-6 - 1.7601372934039951100e-6j)
    assert hermite(-9.5, 100).ae(1.3776300722767084162e-22, abs_eps=0, rel_eps=eps)
    assert hermite(-9.5, -100).ae('1.3106082028470671626e4355')
    assert hermite(-9.5, 100j).ae(-9.7900218581864768430e-23 - 9.7900218581864768430e-23j, abs_eps=0, rel_eps=eps)
    assert hermite(-9.5, -100j).ae(-9.7900218581864768430e-23 + 9.7900218581864768430e-23j, abs_eps=0, rel_eps=eps)
    assert hermite(2+3j, -1-j).ae(851.3677063883687676 - 1496.4373467871007997j)

def test_gegenbauer():
    mp.dps = 15
    assert gegenbauer(1,2,3).ae(12)
    assert gegenbauer(2,3,4).ae(381)
    assert gegenbauer(0,0,0) == 0
    assert gegenbauer(2,-1,3) == 0
    assert gegenbauer(-7, 0.5, 3).ae(8989)
    assert gegenbauer(1, -0.5, 3).ae(-3)
    assert gegenbauer(1, -1.5, 3).ae(-9)
    assert gegenbauer(1, -0.5, 3).ae(-3)
    assert gegenbauer(-0.5, -0.5, 3).ae(-2.6383553159023906245)
    assert gegenbauer(2+3j, 1-j, 3+4j).ae(14.880536623203696780 + 20.022029711598032898j)
    #assert gegenbauer(-2, -0.5, 3).ae(-12)

def test_legenp():
    mp.dps = 15
    assert legenp(2,0,4) == legendre(2,4)
    assert legenp(-2, -1, 0.5).ae(0.43301270189221932338)
    assert legenp(-2, -1, 0.5, type=3).ae(0.43301270189221932338j)
    assert legenp(-2, 1, 0.5).ae(-0.86602540378443864676)
    assert legenp(2+j, 3+4j, -j).ae(134742.98773236786148 + 429782.72924463851745j)
    assert legenp(2+j, 3+4j, -j, type=3).ae(802.59463394152268507 - 251.62481308942906447j)
    assert legenp(2,4,3).ae(0)
    assert legenp(2,4,3,type=3).ae(0)
    assert legenp(2,1,0.5).ae(-1.2990381056766579701)
    assert legenp(2,1,0.5,type=3).ae(1.2990381056766579701j)
    assert legenp(3,2,3).ae(-360)
    assert legenp(3,3,3).ae(240j*2**0.5)
    assert legenp(3,4,3).ae(0)
    assert legenp(0,0.5,2).ae(0.52503756790433198939 - 0.52503756790433198939j)
    assert legenp(-1,-0.5,2).ae(0.60626116232846498110 + 0.60626116232846498110j)
    assert legenp(-2,0.5,2).ae(1.5751127037129959682 - 1.5751127037129959682j)
    assert legenp(-2,0.5,-0.5).ae(-0.85738275810499171286)

def test_legenq():
    mp.dps = 15
    f = legenq
    # Evaluation at poles
    assert isnan(f(3,2,1))
    assert isnan(f(3,2,-1))
    assert isnan(f(3,2,1,type=3))
    assert isnan(f(3,2,-1,type=3))
    # Evaluation at 0
    assert f(0,1,0,type=2).ae(-1)
    assert f(-2,2,0,type=2,zeroprec=200).ae(0)
    assert f(1.5,3,0,type=2).ae(-2.2239343475841951023)
    assert f(0,1,0,type=3).ae(j)
    assert f(-2,2,0,type=3,zeroprec=200).ae(0)
    assert f(1.5,3,0,type=3).ae(2.2239343475841951022*(1-1j))
    # Standard case, degree 0
    assert f(0,0,-1.5).ae(-0.8047189562170501873 + 1.5707963267948966192j)
    assert f(0,0,-0.5).ae(-0.54930614433405484570)
    assert f(0,0,0,zeroprec=200).ae(0)
    assert f(0,0,0.5).ae(0.54930614433405484570)
    assert f(0,0,1.5).ae(0.8047189562170501873 - 1.5707963267948966192j)
    assert f(0,0,-1.5,type=3).ae(-0.80471895621705018730)
    assert f(0,0,-0.5,type=3).ae(-0.5493061443340548457 - 1.5707963267948966192j)
    assert f(0,0,0,type=3).ae(-1.5707963267948966192j)
    assert f(0,0,0.5,type=3).ae(0.5493061443340548457 - 1.5707963267948966192j)
    assert f(0,0,1.5,type=3).ae(0.80471895621705018730)
    # Standard case, degree 1
    assert f(1,0,-1.5).ae(0.2070784343255752810 - 2.3561944901923449288j)
    assert f(1,0,-0.5).ae(-0.72534692783297257715)
    assert f(1,0,0).ae(-1)
    assert f(1,0,0.5).ae(-0.72534692783297257715)
    assert f(1,0,1.5).ae(0.2070784343255752810 - 2.3561944901923449288j)
    # Standard case, degree 2
    assert f(2,0,-1.5).ae(-0.0635669991240192885 + 4.5160394395353277803j)
    assert f(2,0,-0.5).ae(0.81866326804175685571)
    assert f(2,0,0,zeroprec=200).ae(0)
    assert f(2,0,0.5).ae(-0.81866326804175685571)
    assert f(2,0,1.5).ae(0.0635669991240192885 - 4.5160394395353277803j)
    # Misc orders and degrees
    assert f(2,3,1.5,type=2).ae(-5.7243340223994616228j)
    assert f(2,3,1.5,type=3).ae(-5.7243340223994616228)
    assert f(2,3,0.5,type=2).ae(-12.316805742712016310)
    assert f(2,3,0.5,type=3).ae(-12.316805742712016310j)
    assert f(2,3,-1.5,type=2).ae(-5.7243340223994616228j)
    assert f(2,3,-1.5,type=3).ae(5.7243340223994616228)
    assert f(2,3,-0.5,type=2).ae(-12.316805742712016310)
    assert f(2,3,-0.5,type=3).ae(-12.316805742712016310j)
    assert f(2+3j, 3+4j, 0.5, type=3).ae(0.0016119404873235186807 - 0.0005885900510718119836j)
    assert f(2+3j, 3+4j, -1.5, type=3).ae(0.008451400254138808670 + 0.020645193304593235298j)
    assert f(-2.5,1,-1.5).ae(3.9553395527435335749j)
    assert f(-2.5,1,-0.5).ae(1.9290561746445456908)
    assert f(-2.5,1,0).ae(1.2708196271909686299)
    assert f(-2.5,1,0.5).ae(-0.31584812990742202869)
    assert f(-2.5,1,1.5).ae(-3.9553395527435335742 + 0.2993235655044701706j)
    assert f(-2.5,1,-1.5,type=3).ae(0.29932356550447017254j)
    assert f(-2.5,1,-0.5,type=3).ae(-0.3158481299074220287 - 1.9290561746445456908j)
    assert f(-2.5,1,0,type=3).ae(1.2708196271909686292 - 1.2708196271909686299j)
    assert f(-2.5,1,0.5,type=3).ae(1.9290561746445456907 + 0.3158481299074220287j)
    assert f(-2.5,1,1.5,type=3).ae(-0.29932356550447017254)

def test_agm():
    mp.dps = 15
    assert agm(0,0) == 0
    assert agm(0,1) == 0
    assert agm(1,1) == 1
    assert agm(7,7) == 7
    assert agm(j,j) == j
    assert (1/agm(1,sqrt(2))).ae(0.834626841674073186)
    assert agm(1,2).ae(1.4567910310469068692)
    assert agm(1,3).ae(1.8636167832448965424)
    assert agm(1,j).ae(0.599070117367796104+0.599070117367796104j)
    assert agm(2) == agm(1,2)
    assert agm(-3,4).ae(0.63468509766550907+1.3443087080896272j)

def test_gammainc():
    mp.dps = 15
    assert gammainc(2,5).ae(6*exp(-5))
    assert gammainc(2,0,5).ae(1-6*exp(-5))
    assert gammainc(2,3,5).ae(-6*exp(-5)+4*exp(-3))
    assert gammainc(-2.5,-0.5).ae(-0.9453087204829418812-5.3164237738936178621j)
    assert gammainc(0,2,4).ae(0.045121158298212213088)
    assert gammainc(0,3).ae(0.013048381094197037413)
    assert gammainc(0,2+j,1-j).ae(0.00910653685850304839-0.22378752918074432574j)
    assert gammainc(0,1-j).ae(0.00028162445198141833+0.17932453503935894015j)
    assert gammainc(3,4,5,True).ae(0.11345128607046320253)
    assert gammainc(3.5,0,inf).ae(gamma(3.5))
    assert gammainc(-150.5,500).ae('6.9825435345798951153e-627')
    assert gammainc(-150.5,800).ae('4.6885137549474089431e-788')
    assert gammainc(-3.5, -20.5).ae(0.27008820585226911 - 1310.31447140574997636j)
    assert gammainc(-3.5, -200.5).ae(0.27008820585226911 - 5.3264597096208368435e76j) # XXX real part
    assert gammainc(0,0,2) == inf
    assert gammainc(1,b=1).ae(0.6321205588285576784)
    assert gammainc(3,2,2) == 0
    assert gammainc(2,3+j,3-j).ae(-0.28135485191849314194j)
    assert gammainc(4+0j,1).ae(5.8860710587430771455)
    # Regularized upper gamma
    assert isnan(gammainc(0, 0, regularized=True))
    assert gammainc(-1, 0, regularized=True) == inf
    assert gammainc(1, 0, regularized=True) == 1
    assert gammainc(0, 5, regularized=True) == 0
    assert gammainc(0, 2+3j, regularized=True) == 0
    assert gammainc(0, 5000, regularized=True) == 0
    assert gammainc(0, 10**30, regularized=True) == 0
    assert gammainc(-1, 5, regularized=True) == 0
    assert gammainc(-1, 5000, regularized=True) == 0
    assert gammainc(-1, 10**30, regularized=True) == 0
    assert gammainc(-1, -5, regularized=True) == 0
    assert gammainc(-1, -5000, regularized=True) == 0
    assert gammainc(-1, -10**30, regularized=True) == 0
    assert gammainc(-1, 3+4j, regularized=True) == 0
    assert gammainc(1, 5, regularized=True).ae(exp(-5))
    assert gammainc(1, 5000, regularized=True).ae(exp(-5000))
    assert gammainc(1, 10**30, regularized=True).ae(exp(-10**30))
    assert gammainc(1, 3+4j, regularized=True).ae(exp(-3-4j))
    assert gammainc(-1000000,2).ae('1.3669297209397347754e-301037', abs_eps=0, rel_eps=8*eps)
    assert gammainc(-1000000,2,regularized=True) == 0
    assert gammainc(-1000000,3+4j).ae('-1.322575609404222361e-698979 - 4.9274570591854533273e-698978j', abs_eps=0, rel_eps=8*eps)
    assert gammainc(-1000000,3+4j,regularized=True) == 0
    assert gammainc(2+3j, 4+5j, regularized=True).ae(0.085422013530993285774-0.052595379150390078503j)
    assert gammainc(1000j, 1000j, regularized=True).ae(0.49702647628921131761 + 0.00297355675013575341j)
    # Generalized
    assert gammainc(3,4,2) == -gammainc(3,2,4)
    assert gammainc(4, 2, 3).ae(1.2593494302978947396)
    assert gammainc(4, 2, 3, regularized=True).ae(0.20989157171631578993)
    assert gammainc(0, 2, 3).ae(0.035852129613864082155)
    assert gammainc(0, 2, 3, regularized=True) == 0
    assert gammainc(-1, 2, 3).ae(0.015219822548487616132)
    assert gammainc(-1, 2, 3, regularized=True) == 0
    assert gammainc(0, 2, 3).ae(0.035852129613864082155)
    assert gammainc(0, 2, 3, regularized=True) == 0
    # Should use upper gammas
    assert gammainc(5, 10000, 12000).ae('1.1359381951461801687e-4327', abs_eps=0, rel_eps=8*eps)
    # Should use lower gammas
    assert gammainc(10000, 2, 3).ae('8.1244514125995785934e4765')

def test_gammainc_expint_n():
    # These tests are intended to check all cases of the low-level code
    # for upper gamma and expint with small integer index.
    # Need to cover positive/negative arguments; small/large/huge arguments
    # for both positive and negative indices, as well as indices 0 and 1
    # which may be special-cased
    mp.dps = 15
    assert expint(-3,3.5).ae(0.021456366563296693987)
    assert expint(-2,3.5).ae(0.014966633183073309405)
    assert expint(-1,3.5).ae(0.011092916359219041088)
    assert expint(0,3.5).ae(0.0086278238349481430685)
    assert expint(1,3.5).ae(0.0069701398575483929193)
    assert expint(2,3.5).ae(0.0058018939208991255223)
    assert expint(3,3.5).ae(0.0049453773495857807058)
    assert expint(-3,-3.5).ae(-4.6618170604073311319)
    assert expint(-2,-3.5).ae(-5.5996974157555515963)
    assert expint(-1,-3.5).ae(-6.7582555017739415818)
    assert expint(0,-3.5).ae(-9.4615577024835182145)
    assert expint(1,-3.5).ae(-13.925353995152335292 - 3.1415926535897932385j)
    assert expint(2,-3.5).ae(-15.62328702434085977 - 10.995574287564276335j)
    assert expint(3,-3.5).ae(-10.783026313250347722 - 19.242255003237483586j)
    assert expint(-3,350).ae(2.8614825451252838069e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(-2,350).ae(2.8532837224504675901e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(-1,350).ae(2.8451316155828634555e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(0,350).ae(2.8370258275042797989e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(1,350).ae(2.8289659656701459404e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(2,350).ae(2.8209516419468505006e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(3,350).ae(2.8129824725501272171e-155, abs_eps=0, rel_eps=8*eps)
    assert expint(-3,-350).ae(-2.8528796154044839443e+149)
    assert expint(-2,-350).ae(-2.8610072121701264351e+149)
    assert expint(-1,-350).ae(-2.8691813842677537647e+149)
    assert expint(0,-350).ae(-2.8774025343659421709e+149)
    u = expint(1,-350)
    assert u.ae(-2.8856710698020863568e+149)
    assert u.imag.ae(-3.1415926535897932385)
    u = expint(2,-350)
    assert u.ae(-2.8939874026504650534e+149)
    assert u.imag.ae(-1099.5574287564276335)
    u = expint(3,-350)
    assert u.ae(-2.9023519497915044349e+149)
    assert u.imag.ae(-192422.55003237483586)
    assert expint(-3,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(-2,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(-1,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(0,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(1,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(2,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(3,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert expint(-3,-350000000000000000000000).ae('-3.7805306852415755699e+152003068666138139677871')
    assert expint(-2,-350000000000000000000000).ae('-3.7805306852415755699e+152003068666138139677871')
    assert expint(-1,-350000000000000000000000).ae('-3.7805306852415755699e+152003068666138139677871')
    assert expint(0,-350000000000000000000000).ae('-3.7805306852415755699e+152003068666138139677871')
    u = expint(1,-350000000000000000000000)
    assert u.ae('-3.7805306852415755699e+152003068666138139677871')
    assert u.imag.ae(-3.1415926535897932385)
    u = expint(2,-350000000000000000000000)
    assert u.imag.ae(-1.0995574287564276335e+24)
    assert u.ae('-3.7805306852415755699e+152003068666138139677871')
    u = expint(3,-350000000000000000000000)
    assert u.imag.ae(-1.9242255003237483586e+47)
    assert u.ae('-3.7805306852415755699e+152003068666138139677871')
    # Small case; no branch cut
    assert gammainc(-3,3.5).ae(0.00010020262545203707109)
    assert gammainc(-2,3.5).ae(0.00040370427343557393517)
    assert gammainc(-1,3.5).ae(0.0016576839773997501492)
    assert gammainc(0,3.5).ae(0.0069701398575483929193)
    assert gammainc(1,3.5).ae(0.03019738342231850074)
    assert gammainc(2,3.5).ae(0.13588822540043325333)
    assert gammainc(3,3.5).ae(0.64169439772426814072)
    # Small case; with branch cut
    assert gammainc(-3,-3.5).ae(0.03595832954467563286 - 0.52359877559829887308j)
    assert gammainc(-2,-3.5).ae(-0.88024704597962022221 - 1.5707963267948966192j)
    assert gammainc(-1,-3.5).ae(4.4637962926688170771 - 3.1415926535897932385j)
    assert gammainc(0,-3.5).ae(-13.925353995152335292 - 3.1415926535897932385j)
    assert gammainc(1,-3.5).ae(33.115451958692313751)
    assert gammainc(2,-3.5).ae(-82.788629896730784377)
    assert gammainc(3,-3.5).ae(240.08702670051927469)
    # Asymptotic case; no branch cut
    assert gammainc(-3,350).ae(6.5424095113340358813e-163, abs_eps=0, rel_eps=8*eps)
    assert gammainc(-2,350).ae(2.296312222489899769e-160, abs_eps=0, rel_eps=8*eps)
    assert gammainc(-1,350).ae(8.059861834133858573e-158, abs_eps=0, rel_eps=8*eps)
    assert gammainc(0,350).ae(2.8289659656701459404e-155, abs_eps=0, rel_eps=8*eps)
    assert gammainc(1,350).ae(9.9295903962649792963e-153, abs_eps=0, rel_eps=8*eps)
    assert gammainc(2,350).ae(3.485286229089007733e-150, abs_eps=0, rel_eps=8*eps)
    assert gammainc(3,350).ae(1.2233453960006379793e-147, abs_eps=0, rel_eps=8*eps)
    # Asymptotic case; branch cut
    u = gammainc(-3,-350)
    assert u.ae(6.7889565783842895085e+141)
    assert u.imag.ae(-0.52359877559829887308)
    u = gammainc(-2,-350)
    assert u.ae(-2.3692668977889832121e+144)
    assert u.imag.ae(-1.5707963267948966192)
    u = gammainc(-1,-350)
    assert u.ae(8.2685354361441858669e+146)
    assert u.imag.ae(-3.1415926535897932385)
    u = gammainc(0,-350)
    assert u.ae(-2.8856710698020863568e+149)
    assert u.imag.ae(-3.1415926535897932385)
    u = gammainc(1,-350)
    assert u.ae(1.0070908870280797598e+152)
    assert u.imag == 0
    u = gammainc(2,-350)
    assert u.ae(-3.5147471957279983618e+154)
    assert u.imag == 0
    u = gammainc(3,-350)
    assert u.ae(1.2266568422179417091e+157)
    assert u.imag == 0
    # Extreme asymptotic case
    assert gammainc(-3,350000000000000000000000).ae('5.0362468738874738859e-152003068666138139677990', abs_eps=0, rel_eps=8*eps)
    assert gammainc(-2,350000000000000000000000).ae('1.7626864058606158601e-152003068666138139677966', abs_eps=0, rel_eps=8*eps)
    assert gammainc(-1,350000000000000000000000).ae('6.1694024205121555102e-152003068666138139677943', abs_eps=0, rel_eps=8*eps)
    assert gammainc(0,350000000000000000000000).ae('2.1592908471792544286e-152003068666138139677919', abs_eps=0, rel_eps=8*eps)
    assert gammainc(1,350000000000000000000000).ae('7.5575179651273905e-152003068666138139677896', abs_eps=0, rel_eps=8*eps)
    assert gammainc(2,350000000000000000000000).ae('2.645131287794586675e-152003068666138139677872', abs_eps=0, rel_eps=8*eps)
    assert gammainc(3,350000000000000000000000).ae('9.2579595072810533625e-152003068666138139677849', abs_eps=0, rel_eps=8*eps)
    u = gammainc(-3,-350000000000000000000000)
    assert u.ae('8.8175642804468234866e+152003068666138139677800')
    assert u.imag.ae(-0.52359877559829887308)
    u = gammainc(-2,-350000000000000000000000)
    assert u.ae('-3.0861474981563882203e+152003068666138139677824')
    assert u.imag.ae(-1.5707963267948966192)
    u = gammainc(-1,-350000000000000000000000)
    assert u.ae('1.0801516243547358771e+152003068666138139677848')
    assert u.imag.ae(-3.1415926535897932385)
    u = gammainc(0,-350000000000000000000000)
    assert u.ae('-3.7805306852415755699e+152003068666138139677871')
    assert u.imag.ae(-3.1415926535897932385)
    assert gammainc(1,-350000000000000000000000).ae('1.3231857398345514495e+152003068666138139677895')
    assert gammainc(2,-350000000000000000000000).ae('-4.6311500894209300731e+152003068666138139677918')
    assert gammainc(3,-350000000000000000000000).ae('1.6209025312973255256e+152003068666138139677942')

def test_incomplete_beta():
    mp.dps = 15
    assert betainc(-2,-3,0.5,0.75).ae(63.4305673311255413583969)
    assert betainc(4.5,0.5+2j,2.5,6).ae(0.2628801146130621387903065 + 0.5162565234467020592855378j)
    assert betainc(4,5,0,6).ae(90747.77142857142857142857)

def test_erf():
    mp.dps = 15
    assert erf(0) == 0
    assert erf(1).ae(0.84270079294971486934)
    assert erf(3+4j).ae(-120.186991395079444098 - 27.750337293623902498j)
    assert erf(-4-3j).ae(-0.99991066178539168236 + 0.00004972026054496604j)
    assert erf(pi).ae(0.99999112385363235839)
    assert erf(1j).ae(1.6504257587975428760j)
    assert erf(-1j).ae(-1.6504257587975428760j)
    assert isinstance(erf(1), mpf)
    assert isinstance(erf(-1), mpf)
    assert isinstance(erf(0), mpf)
    assert isinstance(erf(0j), mpc)
    assert erf(inf) == 1
    assert erf(-inf) == -1
    assert erfi(0) == 0
    assert erfi(1/pi).ae(0.371682698493894314)
    assert erfi(inf) == inf
    assert erfi(-inf) == -inf
    assert erf(1+0j) == erf(1)
    assert erfc(1+0j) == erfc(1)
    assert erf(0.2+0.5j).ae(1 - erfc(0.2+0.5j))
    assert erfc(0) == 1
    assert erfc(1).ae(1-erf(1))
    assert erfc(-1).ae(1-erf(-1))
    assert erfc(1/pi).ae(1-erf(1/pi))
    assert erfc(-10) == 2
    assert erfc(-1000000) == 2
    assert erfc(-inf) == 2
    assert erfc(inf) == 0
    assert isnan(erfc(nan))
    assert (erfc(10**4)*mpf(10)**43429453).ae('3.63998738656420')
    assert erf(8+9j).ae(-1072004.2525062051158 + 364149.91954310255423j)
    assert erfc(8+9j).ae(1072005.2525062051158 - 364149.91954310255423j)
    assert erfc(-8-9j).ae(-1072003.2525062051158 + 364149.91954310255423j)
    mp.dps = 50
    # This one does not use the asymptotic series
    assert (erfc(10)*10**45).ae('2.0884875837625447570007862949577886115608181193212')
    # This one does
    assert (erfc(50)*10**1088).ae('2.0709207788416560484484478751657887929322509209954')
    mp.dps = 15
    assert str(erfc(10**50)) == '3.66744826532555e-4342944819032518276511289189166050822943970058036665661144537831658646492088707747292249493384317534'
    assert erfinv(0) == 0
    assert erfinv(0.5).ae(0.47693627620446987338)
    assert erfinv(-0.5).ae(-0.47693627620446987338)
    assert erfinv(1) == inf
    assert erfinv(-1) == -inf
    assert erf(erfinv(0.95)).ae(0.95)
    assert erf(erfinv(0.999999999995)).ae(0.999999999995)
    assert erf(erfinv(-0.999999999995)).ae(-0.999999999995)
    mp.dps = 50
    assert erf(erfinv('0.99999999999999999999999999999995')).ae('0.99999999999999999999999999999995')
    assert erf(erfinv('0.999999999999999999999999999999995')).ae('0.999999999999999999999999999999995')
    assert erf(erfinv('-0.999999999999999999999999999999995')).ae('-0.999999999999999999999999999999995')
    mp.dps = 15
    # Complex asymptotic expansions
    v = erfc(50j)
    assert v.real == 1
    assert v.imag.ae('-6.1481820666053078736e+1083')
    assert erfc(-100+5j).ae(2)
    assert (erfc(100+5j)*10**4335).ae(2.3973567853824133572 - 3.9339259530609420597j)
    assert erfc(100+100j).ae(0.00065234366376857698698 - 0.0039357263629214118437j)

def test_pdf():
    mp.dps = 15
    assert npdf(-inf) == 0
    assert npdf(inf) == 0
    assert npdf(5,0,2).ae(npdf(5+4,4,2))
    assert quadts(lambda x: npdf(x,-0.5,0.8), [-inf, inf]) == 1
    assert ncdf(0) == 0.5
    assert ncdf(3,3) == 0.5
    assert ncdf(-inf) == 0
    assert ncdf(inf) == 1
    assert ncdf(10) == 1
    # Verify that this is computed accurately
    assert (ncdf(-10)*10**24).ae(7.619853024160526)

def test_lambertw():
    mp.dps = 15
    assert lambertw(0) == 0
    assert lambertw(0+0j) == 0
    assert lambertw(inf) == inf
    assert isnan(lambertw(nan))
    assert lambertw(inf,1).real == inf
    assert lambertw(inf,1).imag.ae(2*pi)
    assert lambertw(-inf,1).real == inf
    assert lambertw(-inf,1).imag.ae(3*pi)
    assert lambertw(0,-1) == -inf
    assert lambertw(0,1) == -inf
    assert lambertw(0,3) == -inf
    assert lambertw(e).ae(1)
    assert lambertw(1).ae(0.567143290409783873)
    assert lambertw(-pi/2).ae(j*pi/2)
    assert lambertw(-log(2)/2).ae(-log(2))
    assert lambertw(0.25).ae(0.203888354702240164)
    assert lambertw(-0.25).ae(-0.357402956181388903)
    assert lambertw(-1./10000,0).ae(-0.000100010001500266719)
    assert lambertw(-0.25,-1).ae(-2.15329236411034965)
    assert lambertw(0.25,-1).ae(-3.00899800997004620-4.07652978899159763j)
    assert lambertw(-0.25,-1).ae(-2.15329236411034965)
    assert lambertw(0.25,1).ae(-3.00899800997004620+4.07652978899159763j)
    assert lambertw(-0.25,1).ae(-3.48973228422959210+7.41405453009603664j)
    assert lambertw(-4).ae(0.67881197132094523+1.91195078174339937j)
    assert lambertw(-4,1).ae(-0.66743107129800988+7.76827456802783084j)
    assert lambertw(-4,-1).ae(0.67881197132094523-1.91195078174339937j)
    assert lambertw(1000).ae(5.24960285240159623)
    assert lambertw(1000,1).ae(4.91492239981054535+5.44652615979447070j)
    assert lambertw(1000,-1).ae(4.91492239981054535-5.44652615979447070j)
    assert lambertw(1000,5).ae(3.5010625305312892+29.9614548941181328j)
    assert lambertw(3+4j).ae(1.281561806123775878+0.533095222020971071j)
    assert lambertw(-0.4+0.4j).ae(-0.10396515323290657+0.61899273315171632j)
    assert lambertw(3+4j,1).ae(-0.11691092896595324+5.61888039871282334j)
    assert lambertw(3+4j,-1).ae(0.25856740686699742-3.85211668616143559j)
    assert lambertw(-0.5,-1).ae(-0.794023632344689368-0.770111750510379110j)
    assert lambertw(-1./10000,1).ae(-11.82350837248724344+6.80546081842002101j)
    assert lambertw(-1./10000,-1).ae(-11.6671145325663544)
    assert lambertw(-1./10000,-2).ae(-11.82350837248724344-6.80546081842002101j)
    assert lambertw(-1./100000,4).ae(-14.9186890769540539+26.1856750178782046j)
    assert lambertw(-1./100000,5).ae(-15.0931437726379218666+32.5525721210262290086j)
    assert lambertw((2+j)/10).ae(0.173704503762911669+0.071781336752835511j)
    assert lambertw((2+j)/10,1).ae(-3.21746028349820063+4.56175438896292539j)
    assert lambertw((2+j)/10,-1).ae(-3.03781405002993088-3.53946629633505737j)
    assert lambertw((2+j)/10,4).ae(-4.6878509692773249+23.8313630697683291j)
    assert lambertw(-(2+j)/10).ae(-0.226933772515757933-0.164986470020154580j)
    assert lambertw(-(2+j)/10,1).ae(-2.43569517046110001+0.76974067544756289j)
    assert lambertw(-(2+j)/10,-1).ae(-3.54858738151989450-6.91627921869943589j)
    assert lambertw(-(2+j)/10,4).ae(-4.5500846928118151+20.6672982215434637j)
    mp.dps = 50
    assert lambertw(pi).ae('1.073658194796149172092178407024821347547745350410314531')
    mp.dps = 15
    # Former bug in generated branch
    assert lambertw(-0.5+0.002j).ae(-0.78917138132659918344 + 0.76743539379990327749j)
    assert lambertw(-0.5-0.002j).ae(-0.78917138132659918344 - 0.76743539379990327749j)
    assert lambertw(-0.448+0.4j).ae(-0.11855133765652382241 + 0.66570534313583423116j)
    assert lambertw(-0.448-0.4j).ae(-0.11855133765652382241 - 0.66570534313583423116j)
    assert lambertw(-0.65475+0.0001j).ae(-0.61053421111385310898+1.0396534993944097723803j)
    # Huge branch index
    w = lambertw(1,10**20)
    assert w.real.ae(-47.889578926290259164)
    assert w.imag.ae(6.2831853071795864769e+20)

def test_lambertw_hard():
    def check(x,y):
        y = convert(y)
        type_ok = True
        if isinstance(y, mpf):
            type_ok = isinstance(x, mpf)
        real_ok = abs(x.real-y.real) <= abs(y.real)*8*eps
        imag_ok = abs(x.imag-y.imag) <= abs(y.imag)*8*eps
        #print x, y, abs(x.real-y.real), abs(x.imag-y.imag)
        return real_ok and imag_ok
    # Evaluation near 0
    mp.dps = 15
    assert check(lambertw(1e-10), 9.999999999000000000e-11)
    assert check(lambertw(-1e-10), -1.000000000100000000e-10)
    assert check(lambertw(1e-10j), 9.999999999999999999733e-21 + 9.99999999999999999985e-11j)
    assert check(lambertw(-1e-10j), 9.999999999999999999733e-21 - 9.99999999999999999985e-11j)
    assert check(lambertw(1e-10,1), -26.303186778379041559 + 3.265093911703828397j)
    assert check(lambertw(-1e-10,1), -26.326236166739163892 + 6.526183280686333315j)
    assert check(lambertw(1e-10j,1), -26.312931726911421551 + 4.896366881798013421j)
    assert check(lambertw(-1e-10j,1), -26.297238779529035066 + 1.632807161345576513j)
    assert check(lambertw(1e-10,-1), -26.303186778379041559 - 3.265093911703828397j)
    assert check(lambertw(-1e-10,-1), -26.295238819246925694)
    assert check(lambertw(1e-10j,-1), -26.297238779529035028 - 1.6328071613455765135j)
    assert check(lambertw(-1e-10j,-1), -26.312931726911421551 - 4.896366881798013421j)
    # Test evaluation very close to the branch point -1/e
    # on the -1, 0, and 1 branches
    add = lambda x, y: fadd(x,y,exact=True)
    sub = lambda x, y: fsub(x,y,exact=True)
    addj = lambda x, y: fadd(x,fmul(y,1j,exact=True),exact=True)
    subj = lambda x, y: fadd(x,fmul(y,-1j,exact=True),exact=True)
    mp.dps = 1500
    a = -1/e + 10*eps
    d3 = mpf('1e-3')
    d10 = mpf('1e-10')
    d20 = mpf('1e-20')
    d40 = mpf('1e-40')
    d80 = mpf('1e-80')
    d300 = mpf('1e-300')
    d1000 = mpf('1e-1000')
    mp.dps = 15
    # ---- Branch 0 ----
    # -1/e + eps
    assert check(lambertw(add(a,d3)), -0.92802015005456704876)
    assert check(lambertw(add(a,d10)), -0.99997668374140088071)
    assert check(lambertw(add(a,d20)), -0.99999999976683560186)
    assert lambertw(add(a,d40)) == -1
    assert lambertw(add(a,d80)) == -1
    assert lambertw(add(a,d300)) == -1
    assert lambertw(add(a,d1000)) == -1
    # -1/e - eps
    assert check(lambertw(sub(a,d3)), -0.99819016149860989001+0.07367191188934638577j)
    assert check(lambertw(sub(a,d10)), -0.9999999998187812114595992+0.0000233164398140346109194j)
    assert check(lambertw(sub(a,d20)), -0.99999999999999999998187+2.331643981597124203344e-10j)
    assert check(lambertw(sub(a,d40)), -1.0+2.33164398159712420336e-20j)
    assert check(lambertw(sub(a,d80)), -1.0+2.33164398159712420336e-40j)
    assert check(lambertw(sub(a,d300)), -1.0+2.33164398159712420336e-150j)
    assert check(lambertw(sub(a,d1000)), mpc(-1,'2.33164398159712420336e-500'))
    # -1/e + eps*j
    assert check(lambertw(addj(a,d3)), -0.94790387486938526634+0.05036819639190132490j)
    assert check(lambertw(addj(a,d10)), -0.9999835127872943680999899+0.0000164870314895821225256j)
    assert check(lambertw(addj(a,d20)), -0.999999999835127872929987+1.64872127051890935830e-10j)
    assert check(lambertw(addj(a,d40)), -0.9999999999999999999835+1.6487212707001281468305e-20j)
    assert check(lambertw(addj(a,d80)), -1.0 + 1.64872127070012814684865e-40j)
    assert check(lambertw(addj(a,d300)), -1.0 + 1.64872127070012814684865e-150j)
    assert check(lambertw(addj(a,d1000)), mpc(-1.0,'1.64872127070012814684865e-500'))
    # -1/e - eps*j
    assert check(lambertw(subj(a,d3)), -0.94790387486938526634-0.05036819639190132490j)
    assert check(lambertw(subj(a,d10)), -0.9999835127872943680999899-0.0000164870314895821225256j)
    assert check(lambertw(subj(a,d20)), -0.999999999835127872929987-1.64872127051890935830e-10j)
    assert check(lambertw(subj(a,d40)), -0.9999999999999999999835-1.6487212707001281468305e-20j)
    assert check(lambertw(subj(a,d80)), -1.0 - 1.64872127070012814684865e-40j)
    assert check(lambertw(subj(a,d300)), -1.0 - 1.64872127070012814684865e-150j)
    assert check(lambertw(subj(a,d1000)), mpc(-1.0,'-1.64872127070012814684865e-500'))
    # ---- Branch 1 ----
    assert check(lambertw(addj(a,d3),1), -3.088501303219933378005990 + 7.458676867597474813950098j)
    assert check(lambertw(addj(a,d80),1), -3.088843015613043855957087 + 7.461489285654254556906117j)
    assert check(lambertw(addj(a,d300),1), -3.088843015613043855957087 + 7.461489285654254556906117j)
    assert check(lambertw(addj(a,d1000),1), -3.088843015613043855957087 + 7.461489285654254556906117j)
    assert check(lambertw(subj(a,d3),1), -1.0520914180450129534365906 + 0.0539925638125450525673175j)
    assert check(lambertw(subj(a,d10),1), -1.0000164872127056318529390 + 0.000016487393927159250398333077j)
    assert check(lambertw(subj(a,d20),1), -1.0000000001648721270700128 + 1.64872127088134693542628e-10j)
    assert check(lambertw(subj(a,d40),1), -1.000000000000000000016487 + 1.64872127070012814686677e-20j)
    assert check(lambertw(subj(a,d80),1), -1.0 + 1.64872127070012814684865e-40j)
    assert check(lambertw(subj(a,d300),1), -1.0 + 1.64872127070012814684865e-150j)
    assert check(lambertw(subj(a,d1000),1), mpc(-1.0, '1.64872127070012814684865e-500'))
    # ---- Branch -1 ----
    # -1/e + eps
    assert check(lambertw(add(a,d3),-1), -1.075608941186624989414945)
    assert check(lambertw(add(a,d10),-1), -1.000023316621036696460620)
    assert check(lambertw(add(a,d20),-1), -1.000000000233164398177834)
    assert lambertw(add(a,d40),-1) == -1
    assert lambertw(add(a,d80),-1) == -1
    assert lambertw(add(a,d300),-1) == -1
    assert lambertw(add(a,d1000),-1) == -1
    # -1/e - eps
    assert check(lambertw(sub(a,d3),-1), -0.99819016149860989001-0.07367191188934638577j)
    assert check(lambertw(sub(a,d10),-1), -0.9999999998187812114595992-0.0000233164398140346109194j)
    assert check(lambertw(sub(a,d20),-1), -0.99999999999999999998187-2.331643981597124203344e-10j)
    assert check(lambertw(sub(a,d40),-1), -1.0-2.33164398159712420336e-20j)
    assert check(lambertw(sub(a,d80),-1), -1.0-2.33164398159712420336e-40j)
    assert check(lambertw(sub(a,d300),-1), -1.0-2.33164398159712420336e-150j)
    assert check(lambertw(sub(a,d1000),-1), mpc(-1,'-2.33164398159712420336e-500'))
    # -1/e + eps*j
    assert check(lambertw(addj(a,d3),-1), -1.0520914180450129534365906 - 0.0539925638125450525673175j)
    assert check(lambertw(addj(a,d10),-1), -1.0000164872127056318529390 - 0.0000164873939271592503983j)
    assert check(lambertw(addj(a,d20),-1), -1.0000000001648721270700 - 1.64872127088134693542628e-10j)
    assert check(lambertw(addj(a,d40),-1), -1.00000000000000000001648 - 1.6487212707001281468667726e-20j)
    assert check(lambertw(addj(a,d80),-1), -1.0 - 1.64872127070012814684865e-40j)
    assert check(lambertw(addj(a,d300),-1), -1.0 - 1.64872127070012814684865e-150j)
    assert check(lambertw(addj(a,d1000),-1), mpc(-1.0,'-1.64872127070012814684865e-500'))
    # -1/e - eps*j
    assert check(lambertw(subj(a,d3),-1), -3.088501303219933378005990-7.458676867597474813950098j)
    assert check(lambertw(subj(a,d10),-1), -3.088843015579260686911033-7.461489285372968780020716j)
    assert check(lambertw(subj(a,d20),-1), -3.088843015613043855953708-7.461489285654254556877988j)
    assert check(lambertw(subj(a,d40),-1), -3.088843015613043855957087-7.461489285654254556906117j)
    assert check(lambertw(subj(a,d80),-1), -3.088843015613043855957087 - 7.461489285654254556906117j)
    assert check(lambertw(subj(a,d300),-1), -3.088843015613043855957087 - 7.461489285654254556906117j)
    assert check(lambertw(subj(a,d1000),-1), -3.088843015613043855957087 - 7.461489285654254556906117j)
    # One more case, testing higher precision
    mp.dps = 500
    x = -1/e + mpf('1e-13')
    ans = "-0.99999926266961377166355784455394913638782494543377383"\
    "744978844374498153493943725364881490261187530235150668593869563"\
    "168276697689459394902153960200361935311512317183678882"
    mp.dps = 15
    assert lambertw(x).ae(ans)
    mp.dps = 50
    assert lambertw(x).ae(ans)
    mp.dps = 150
    assert lambertw(x).ae(ans)

def test_meijerg():
    mp.dps = 15
    assert meijerg([[2,3],[1]],[[0.5,2],[3,4]], 2.5).ae(4.2181028074787439386)
    assert meijerg([[],[1+j]],[[1],[1]], 3+4j).ae(271.46290321152464592 - 703.03330399954820169j)
    assert meijerg([[0.25],[1]],[[0.5],[2]],0) == 0
    assert meijerg([[0],[]],[[0,0,'1/3','2/3'], []], '2/27').ae(2.2019391389653314120)
    # Verify 1/z series being used
    assert meijerg([[-3],[-0.5]], [[-1],[-2.5]], -0.5).ae(-1.338096165935754898687431)
    assert meijerg([[1-(-1)],[1-(-2.5)]], [[1-(-3)],[1-(-0.5)]], -2.0).ae(-1.338096165935754898687431)
    assert meijerg([[-3],[-0.5]], [[-1],[-2.5]], -1).ae(-(pi+4)/(4*pi))
    a = 2.5
    b = 1.25
    for z in [mpf(0.25), mpf(2)]:
        x1 = hyp1f1(a,b,z)
        x2 = gamma(b)/gamma(a)*meijerg([[1-a],[]],[[0],[1-b]],-z)
        x3 = gamma(b)/gamma(a)*meijerg([[1-0],[1-(1-b)]],[[1-(1-a)],[]],-1/z)
        assert x1.ae(x2)
        assert x1.ae(x3)

def test_appellf1():
    mp.dps = 15
    assert appellf1(2,-2,1,1,2,3).ae(-1.75)
    assert appellf1(2,1,-2,1,2,3).ae(-8)
    assert appellf1(2,1,-2,1,0.5,0.25).ae(1.5)
    assert appellf1(-2,1,3,2,3,3).ae(19)
    assert appellf1(1,2,3,4,0.5,0.125).ae( 1.53843285792549786518)

def test_coulomb():
    # Note: most tests are doctests
    # Test for a bug:
    mp.dps = 15
    assert coulombg(mpc(-5,0),2,3).ae(20.087729487721430394)

def test_hyper_param_accuracy():
    mp.dps = 15
    As = [n+1e-10 for n in range(-5,-1)]
    Bs = [n+1e-10 for n in range(-12,-5)]
    assert hyper(As,Bs,10).ae(-381757055858.652671927)
    assert legenp(0.5, 100, 0.25).ae(-2.4124576567211311755e+144)
    assert (hyp1f1(1000,1,-100)*10**24).ae(5.2589445437370169113)
    assert (hyp2f1(10, -900, 10.5, 0.99)*10**24).ae(1.9185370579660768203)
    assert (hyp2f1(1000,1.5,-3.5,-1.5)*10**385).ae(-2.7367529051334000764)
    assert hyp2f1(-5, 10, 3, 0.5, zeroprec=500) == 0
    assert (hyp1f1(-10000, 1000, 100)*10**424).ae(-3.1046080515824859974)
    assert (hyp2f1(1000,1.5,-3.5,-0.75,maxterms=100000)*10**231).ae(-4.0534790813913998643)
    assert legenp(2, 3, 0.25) == 0
    try:
        hypercomb(lambda a: [([],[],[],[],[a],[-a],0.5)], [3])
        assert 0
    except ValueError:
        pass
    assert hypercomb(lambda a: [([],[],[],[],[a],[-a],0.5)], [3], infprec=200) == inf
    assert meijerg([[],[]],[[0,0,0,0],[]],0.1).ae(1.5680822343832351418)
    assert (besselk(400,400)*10**94).ae(1.4387057277018550583)
    mp.dps = 5
    (hyp1f1(-5000.5, 1500, 100)*10**185).ae(8.5185229673381935522)
    (hyp1f1(-5000, 1500, 100)*10**185).ae(9.1501213424563944311)
    mp.dps = 15
    (hyp1f1(-5000.5, 1500, 100)*10**185).ae(8.5185229673381935522)
    (hyp1f1(-5000, 1500, 100)*10**185).ae(9.1501213424563944311)
    assert hyp0f1(fadd(-20,'1e-100',exact=True), 0.25).ae(1.85014429040102783e+49)
    assert hyp0f1((-20*10**100+1, 10**100), 0.25).ae(1.85014429040102783e+49)

def test_hypercomb_zero_pow():
    # check that 0^0 = 1
    assert hypercomb(lambda a: (([0],[a],[],[],[],[],0),), [0]) == 1
    assert meijerg([[-1.5],[]],[[0],[-0.75]],0).ae(1.4464090846320771425)

def test_spherharm():
    mp.dps = 15
    t = 0.5; r = 0.25
    assert spherharm(0,0,t,r).ae(0.28209479177387814347)
    assert spherharm(1,-1,t,r).ae(0.16048941205971996369 - 0.04097967481096344271j)
    assert spherharm(1,0,t,r).ae(0.42878904414183579379)
    assert spherharm(1,1,t,r).ae(-0.16048941205971996369 - 0.04097967481096344271j)
    assert spherharm(2,-2,t,r).ae(0.077915886919031181734 - 0.042565643022253962264j)
    assert spherharm(2,-1,t,r).ae(0.31493387233497459884 - 0.08041582001959297689j)
    assert spherharm(2,0,t,r).ae(0.41330596756220761898)
    assert spherharm(2,1,t,r).ae(-0.31493387233497459884 - 0.08041582001959297689j)
    assert spherharm(2,2,t,r).ae(0.077915886919031181734 + 0.042565643022253962264j)
    assert spherharm(3,-3,t,r).ae(0.033640236589690881646 - 0.031339125318637082197j)
    assert spherharm(3,-2,t,r).ae(0.18091018743101461963 - 0.09883168583167010241j)
    assert spherharm(3,-1,t,r).ae(0.42796713930907320351 - 0.10927795157064962317j)
    assert spherharm(3,0,t,r).ae(0.27861659336351639787)
    assert spherharm(3,1,t,r).ae(-0.42796713930907320351 - 0.10927795157064962317j)
    assert spherharm(3,2,t,r).ae(0.18091018743101461963 + 0.09883168583167010241j)
    assert spherharm(3,3,t,r).ae(-0.033640236589690881646 - 0.031339125318637082197j)
    assert spherharm(0,-1,t,r) == 0
    assert spherharm(0,-2,t,r) == 0
    assert spherharm(0,1,t,r) == 0
    assert spherharm(0,2,t,r) == 0
    assert spherharm(1,2,t,r) == 0
    assert spherharm(1,3,t,r) == 0
    assert spherharm(1,-2,t,r) == 0
    assert spherharm(1,-3,t,r) == 0
    assert spherharm(2,3,t,r) == 0
    assert spherharm(2,4,t,r) == 0
    assert spherharm(2,-3,t,r) == 0
    assert spherharm(2,-4,t,r) == 0
    assert spherharm(3,4.5,0.5,0.25).ae(-22.831053442240790148 + 10.910526059510013757j)
    assert spherharm(2+3j, 1-j, 1+j, 3+4j).ae(-2.6582752037810116935 - 1.0909214905642160211j)
    assert spherharm(-6,2.5,t,r).ae(0.39383644983851448178 + 0.28414687085358299021j)
    assert spherharm(-3.5, 3, 0.5, 0.25).ae(0.014516852987544698924 - 0.015582769591477628495j)
    assert spherharm(-3, 3, 0.5, 0.25) == 0
    assert spherharm(-6, 3, 0.5, 0.25).ae(-0.16544349818782275459 - 0.15412657723253924562j)
    assert spherharm(-6, 1.5, 0.5, 0.25).ae(0.032208193499767402477 + 0.012678000924063664921j)
    assert spherharm(3,0,0,1).ae(0.74635266518023078283)
    assert spherharm(3,-2,0,1) == 0
    assert spherharm(3,-2,1,1).ae(-0.16270707338254028971 - 0.35552144137546777097j)

def test_qfunctions():
    mp.dps = 15
    assert qp(2,3,100).ae('2.7291482267247332183e2391')

def test_issue_3298():
    mp.prec = 150
    x = ldexp(2476979795053773,-52)
    assert betainc(206, 385, 0, 0.55, 1).ae('0.99999999999999999999996570910644857895771110649954')
    mp.dps = 15
    try:
        u = hyp2f1(-5,5,0.5,0.5)
        raise AssertionError("hyp2f1(-5,5,0.5,0.5) (failed zero detection)")
    except (mp.NoConvergence, ValueError):
        pass
