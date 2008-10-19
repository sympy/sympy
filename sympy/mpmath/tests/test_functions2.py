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
    assert jv(5,1).ae(0.000249757730211234431)
    assert jv(5,pi).ae(0.0521411843671184747)
    assert jv(5,1000).ae(0.00502540694523318607)
    assert jv(5,-25).ae(0.0660079953984229934)
    assert jv(-3,2).ae(-0.128943249474402051)
    assert jv(-4,2).ae(0.0339957198075684341)
    assert jv(3,3+2j).ae(0.424718794929639595942 + 0.625665327745785804812j)
    assert jv(0.25,4).ae(-0.374760630804249715)
    assert jv(1+2j,3+4j).ae(0.319247428741872131 - 0.669557748880365678j)

def test_hyper_misc():
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
    assert ellipk(-1).ae('1.3110287771460599052')
    assert ellipk(-2).ae('1.1714200841467698589')
    assert isinstance(ellipk(-2), mpf)
    assert ellipk(-50).ae('0.47103424540873331679')
    assert ellipk(fraction(99999,100000)).ae('7.1427724505817781901')
    assert ellipk(fraction(100001,100000)).ae(mpc('7.1427417367963090109', '-1.5707923998261688019'))
    assert ellipk(2).ae(mpc('1.3110287771460599052', '-1.3110287771460599052'))
    assert ellipk(50).ae(mpc('0.22326753950210985451', '-0.47434723226254522087'))
    assert ellipk(3+4j).ae(mpc('0.91119556380496500866', '0.63133428324134524388'))
    assert ellipk(3-4j).ae(mpc('0.91119556380496500866', '-0.63133428324134524388'))
    assert ellipk(-3+4j).ae(mpc('0.95357894880405122483', '0.23093044503746114444'))
    assert ellipk(-3-4j).ae(mpc('0.95357894880405122483', '-0.23093044503746114444'))
    assert isnan(ellipk(nan))
    assert ellipk(inf) == 0
    assert ellipk(-inf) == 0
    assert ellipe(0).ae(pi/2)
    assert ellipe(0.5).ae(pi**(mpf(3)/2)/gamma(0.25)**2 +gamma(0.25)**2/(8*sqrt(pi)))
    assert ellipe(1) == 1

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
    assert ei(inf) == inf
    assert ei(-inf) == -0.0
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

def test_hyper_0f1():
    v = 8.63911136507950465
    assert hyper([],[(1,3)],1.5).ae(v)
    assert hyper([],[1/3.],1.5).ae(v)
    assert hyp0f1(1/3.,1.5).ae(v)
    assert hyp0f1((1,3),1.5).ae(v)

def test_hyper_1f1():
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

def test_hyper_2f1():
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

def test_orthpoly():
    mp.dps = 15
    assert jacobi(-4,2,3,0.7).ae(22800./4913)
    assert jacobi(3,2,4,5.5) == 4133.125
    assert jacobi(1.5,5/6.,4,0).ae(-1.0851951434075508417)
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

def test_incomplete_gamma():
    mp.dps = 15
    assert upper_gamma(-2.5,-0.5).ae(-0.9453087204829418812-5.3164237738936178621j)

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
    mp.dps = 50
    # This one does not use the asymptotic series
    assert (erfc(10)*10**45).ae('2.0884875837625447570007862949577886115608181193212')
    # This one does
    assert (erfc(50)*10**1088).ae('2.0709207788416560484484478751657887929322509209954')
    mp.dps = 15
    assert str(erfc(10**50)) == '3.66744826532555e-4342944819032518276511289189166050822943970058036665661144537831658646492088707747292249493384317534'

def test_pdf():
    mp.dps = 15
    assert npdf(-inf) == 0
    assert npdf(inf) == 0
    assert npdf(5,0,2).ae(5+4,4,2)
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
    assert isnan(lambertw(-inf))
    assert isnan(lambertw(nan))
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
