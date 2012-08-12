"""Tests for Groebner bases. """

from sympy.polys.distributedpolys import (
    sdp_from_dict,
)

from sympy.polys.groebnertools import (
    sdp_groebner, sig, sig_key, sig_cmp,
    lbp, lbp_cmp, lbp_key, critical_pair,
    cp_cmp, cp_key, is_rewritable_or_comparable,
    Sign, Polyn, Num, s_poly, f5_reduce,
    _basis, _representing_matrices,
    matrix_fglm,
)

from sympy.polys.monomialtools import (
    lex, grlex, grevlex,
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, DomainError,
)

from sympy.polys.domains import ZZ, QQ

from sympy import S, Symbol, symbols, groebner
from sympy.utilities.pytest import raises, skip, XFAIL
from sympy.polys import polyconfig as config

def helper_test_sdp_groebner():
    f = sdp_from_dict({(1,2): QQ(2,), (2,0): QQ(1)}, lex)
    g = sdp_from_dict({(0,3): QQ(2), (1,1): QQ(1), (0,0): QQ(-1)},  lex)

    a = sdp_from_dict({(1,0): QQ(1,1)}, lex)
    b = sdp_from_dict({(0,3): QQ(1,1), (0,0): QQ(-1,2)}, lex)

    assert sdp_groebner((f, g), 1, lex, QQ) == [a, b]

    f = sdp_from_dict({(2,1): QQ(2,), (0,2): QQ(1)}, lex)
    g = sdp_from_dict({(3,0): QQ(2), (1,1): QQ(1), (0,0): QQ(-1)},  lex)

    a = sdp_from_dict({(0,1): QQ(1,1)}, lex)
    b = sdp_from_dict({(3,0): QQ(1,1), (0,0): QQ(-1,2)}, lex)

    assert sdp_groebner((f, g), 1, lex, QQ) == [b, a]

    f = sdp_from_dict({(0,0,2): QQ(-1), (1,0,0): QQ(1)}, lex)
    g = sdp_from_dict({(0,0,3): QQ(-1), (0,1,0): QQ(1)}, lex)

    assert sdp_groebner((f, g), 1, lex, QQ) == [f, g]

    f = sdp_from_dict({(3,0): QQ(1), (1,1): QQ(-2)}, grlex)
    g = sdp_from_dict({(2,1): QQ(1), (0,2): QQ(-2), (1,0): QQ(1)}, grlex)

    a = sdp_from_dict({(2,0): QQ(1)}, grlex)
    b = sdp_from_dict({(1,1): QQ(1)}, grlex)
    c = sdp_from_dict({(0,2): QQ(1), (1, 0): QQ(-1,2)}, grlex)

    assert sdp_groebner((f, g), 1, grlex, QQ) == [a, b, c]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,1,0): QQ(1)}, lex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,0,1): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,0): -QQ(1)}, lex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, lex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, lex),
        sdp_from_dict({(0,3,0): QQ(1), (0,0,2): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,1,0): QQ(1)}, grlex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,0,1): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(0,3,0): QQ(1), (0,0,2): -QQ(1)}, grlex),
        sdp_from_dict({(2,0,0): QQ(1), (0,1,0): -QQ(1)}, grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,0,1): QQ(1)}, lex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,1,0): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,0,1): -QQ(1)}, lex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, lex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, lex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,3): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,0,1): QQ(1)}, grlex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,1,0): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(0,0,3): QQ(1), (0,2,0): -QQ(1)}, grlex),
        sdp_from_dict({(2,0,0): QQ(1), (0,0,1): -QQ(1)}, grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (1,0,0): QQ(1)}, lex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (0,0,1): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,2,0): -QQ(1)}, lex),
        sdp_from_dict({(0,3,0): QQ(1), (0,0,1): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (1,0,0): QQ(1)}, grlex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (0,0,1): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,1): -QQ(1)}, grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, grlex),
        sdp_from_dict({(0,2,0): QQ(1), (1,0,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (1,0,0): QQ(1)}, lex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (0,1,0): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,0,2): -QQ(1)}, lex),
        sdp_from_dict({(0,1,0): QQ(1), (0,0,3): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (1,0,0): QQ(1)}, grlex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (0,1,0): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,1): -QQ(1)}, grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, grlex),
        sdp_from_dict({(0,0,2): QQ(1), (1,0,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (0,0,1): QQ(1)}, lex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (1,0,0): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,1,1): -QQ(1)}, lex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,1): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (0,0,1): QQ(1)}, grlex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (1,0,0): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(0,0,3): QQ(1), (2,0,0): -QQ(1)}, grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, grlex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,1): -QQ(1)}, grlex),
        sdp_from_dict({(0,1,1): QQ(1), (1,0,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (0,1,0): QQ(1)}, lex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (1,0,0): QQ(1)}, lex)

    assert sdp_groebner((f, g), 2, lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,0,3): -QQ(1)}, lex),
        sdp_from_dict({(0,1,0): QQ(1), (0,0,2): -QQ(1)}, lex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (0,1,0): QQ(1)}, grlex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (1,0,0): QQ(1)}, grlex)

    assert sdp_groebner((f, g), 2, grlex, QQ) == [
        sdp_from_dict({(0,3,0): QQ(1), (2,0,0): -QQ(1)}, grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, grlex),
        sdp_from_dict({(0,1,1): QQ(1), (1,0,0): -QQ(1)}, grlex),
        sdp_from_dict({(0,0,2): QQ(1), (0,1,0): -QQ(1)}, grlex),
    ]

    f = sdp_from_dict({(2,2): QQ(4), (1,1): QQ(4), (0,0): QQ(1)}, lex)
    g = sdp_from_dict({(2,0): QQ(1), (0,2): QQ(1), (0,0):-QQ(1)}, lex)

    assert sdp_groebner((f, g), 1, lex, QQ) == [
        sdp_from_dict({(1,0): QQ(1,1), (0,7): QQ(-4,1), (0,5): QQ(8,1), (0,3): QQ(-7,1), (0,1): QQ(3,1)}, lex),
        sdp_from_dict({(0,8): QQ(1,1), (0,6): QQ(-2,1), (0,4): QQ(3,2), (0,2): QQ(-1,2), (0,0): QQ(1,16)}, lex),
    ]

    raises(DomainError, lambda: sdp_groebner([], 1, lex, ZZ))

def test_sdp_groebner():
    config.setup('GB_METHOD', 'f5b')
    helper_test_sdp_groebner()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_sdp_groebner()

def helper_test_benchmark_minpoly():
    x, y, z = symbols('x,y,z')

    I = [x**3 + x + 1, y**2 + y + 1, (x + y) * z - (x**2 + y)]

    assert groebner(I, x, y, z, order='lex') == [
        -975 + 2067*x + 6878*z - 11061*z**2 + 6062*z**3 - 1065*z**4 + 155*z**5,
        -308 + 159*y + 1043*z - 1161*z**2 + 523*z**3 - 91*z**4 + 12*z**5,
        13 - 46*z + 89*z**2 - 82*z**3 + 41*z**4 - 7*z**5 + z**6,
    ]

    assert groebner(I, x, y, z, order='lex', field=True) == [
        -S(25)/53 + x + 6878*z/2067 - 3687*z**2/689 + 6062*z**3/2067 - 355*z**4/689 + 155*z**5/2067,
        -S(308)/159 + y + 1043*z/159 - 387*z**2/53 + 523*z**3/159 - 91*z**4/159 + 4*z**5/53,
        13 - 46*z + 89*z**2 - 82*z**3 + 41*z**4 - 7*z**5 + z**6,
    ]

def test_benchmark_minpoly():
    config.setup('GB_METHOD', 'f5b')
    helper_test_benchmark_minpoly()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_benchmark_minpoly()

@XFAIL
def test_benchmark_coloring():
    skip('takes too much time')

    V = range(1, 12+1)
    E = [(1,2),(2,3),(1,4),(1,6),(1,12),(2,5),(2,7),(3,8),(3,10),
         (4,11),(4,9),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),
         (11,12),(5,12),(5,9),(6,10),(7,11),(8,12),(3,4)]

    V = [Symbol('x' + str(i)) for i in V]
    E = [(V[i-1], V[j-1]) for i, j in E]

    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = V

    I3 = [x**3 - 1 for x in V]
    Ig = [x**2 + x*y + y**2 for x, y in E]

    I = I3 + Ig

    assert groebner(I[:-1], V, order='lex') == [
        x1 + x11 + x12,
        x2 - x11,
        x3 - x12,
        x4 - x12,
        x5 + x11 + x12,
        x6 - x11,
        x7 - x12,
        x8 + x11 + x12,
        x9 - x11,
        x10 + x11 + x12,
        x11**2 + x11*x12 + x12**2,
        x12**3 - 1,
    ]

    assert groebner(I, V, order='lex') == [1]

def helper_test_benchmark_katsura_3():
    x0, x1, x2 = symbols('x:3')

    I = [x0 + 2*x1 + 2*x2 - 1,
         x0**2 + 2*x1**2 + 2*x2**2 - x0,
         2*x0*x1 + 2*x1*x2 - x1]

    assert groebner(I, x0, x1, x2, order='lex') == [
        -7 + 7*x0 + 8*x2 + 158*x2**2 - 420*x2**3,
        7*x1 + 3*x2 - 79*x2**2 + 210*x2**3,
        x2 + x2**2 - 40*x2**3 + 84*x2**4,
    ]

    assert groebner(I, x0, x1, x2, order='grlex') == [
        7*x1 + 3*x2 - 79*x2**2 + 210*x2**3,
        -x1 + x2 - 3*x2**2 + 5*x1**2,
        -x1 - 4*x2 + 10*x1*x2 + 12*x2**2,
        -1 + x0 + 2*x1 + 2*x2,
    ]

def test_benchmark_katsura3():
    config.setup('GB_METHOD', 'f5b')
    helper_test_benchmark_katsura_3()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_benchmark_katsura_3()

def helper_test_benchmark_katsura_4():
    x0, x1, x2, x3 = symbols('x:4')

    I = [x0 + 2*x1 + 2*x2 + 2*x3 - 1,
         x0**2 + 2*x1**2 + 2*x2**2 + 2*x3**2 - x0,
         2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1,
         x1**2 + 2*x0*x2 + 2*x1*x3 - x2]

    assert groebner(I, x0, x1, x2, x3, order='lex') == [
        5913075*x0 - 159690237696*x3**7 + 31246269696*x3**6 + 27439610544*x3**5 - 6475723368*x3**4 - 838935856*x3**3 + 275119624*x3**2 + 4884038*x3 - 5913075,
        1971025*x1 - 97197721632*x3**7 + 73975630752*x3**6 - 12121915032*x3**5 - 2760941496*x3**4 + 814792828*x3**3 - 1678512*x3**2 - 9158924*x3,
        5913075*x2 + 371438283744*x3**7 - 237550027104*x3**6 + 22645939824*x3**5 + 11520686172*x3**4 - 2024910556*x3**3 - 132524276*x3**2 + 30947828*x3,
        128304*x3**8 - 93312*x3**7 + 15552*x3**6 + 3144*x3**5 - 1120*x3**4 + 36*x3**3 + 15*x3**2 - x3,
    ]

    assert groebner(I, x0, x1, x2, x3, order='grlex') == [
        393*x1 - 4662*x2**2 + 4462*x2*x3 - 59*x2 + 224532*x3**4 - 91224*x3**3 - 678*x3**2 + 2046*x3,
        -x1 + 196*x2**3 - 21*x2**2 + 60*x2*x3 - 18*x2 - 168*x3**3 + 83*x3**2 - 9*x3,
        -6*x1 + 1134*x2**2*x3 - 189*x2**2 - 466*x2*x3 + 32*x2 - 630*x3**3 + 57*x3**2 + 51*x3,
        33*x1 + 63*x2**2 + 2268*x2*x3**2 - 188*x2*x3 + 34*x2 + 2520*x3**3 - 849*x3**2 + 3*x3,
        7*x1**2 - x1 - 7*x2**2 - 24*x2*x3 + 3*x2 - 15*x3**2 + 5*x3,
        14*x1*x2 - x1 + 14*x2**2 + 18*x2*x3 - 4*x2 + 6*x3**2 - 2*x3,
        14*x1*x3 - x1 + 7*x2**2 + 32*x2*x3 - 4*x2 + 27*x3**2 - 9*x3,
        x0 + 2*x1 + 2*x2 + 2*x3 - 1,
    ]

def test_benchmark_kastura_4():
    config.setup('GB_METHOD', 'f5b')
    helper_test_benchmark_katsura_4()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_benchmark_katsura_4()

def helper_test_benchmark_czichowski():
    x, t = symbols('x t')

    I = [9*x**8 + 36*x**7 - 32*x**6 - 252*x**5 - 78*x**4 + 468*x**3 + 288*x**2 - 108*x + 9, (-72 - 72*t)*x**7 + (-256 - 252*t)*x**6 + (192 + 192*t)*x**5 + (1280 + 1260*t)*x**4 + (312 + 312*t)*x**3 + (-404*t)*x**2 + (-576 - 576*t)*x + 96 + 108*t]

    assert groebner(I, x, t, order='lex') == [
        -160420835591776763325581422211936558925462474417709511019228211783493866564923546661604487873*t**7 - 1406108495478033395547109582678806497509499966197028487131115097902188374051595011248311352864*t**6 - 5241326875850889518164640374668786338033653548841427557880599579174438246266263602956254030352*t**5 - 10758917262823299139373269714910672770004760114329943852726887632013485035262879510837043892416*t**4 - 13119383576444715672578819534846747735372132018341964647712009275306635391456880068261130581248*t**3 - 9491412317016197146080450036267011389660653495578680036574753839055748080962214787557853941760*t**2 - 3767520915562795326943800040277726397326609797172964377014046018280260848046603967211258368000*t + 3725588592068034903797967297424801242396746870413359539263038139343329273586196480000*x - 632314652371226552085897259159210286886724229880266931574701654721512325555116066073245696000,
        610733380717522355121*t**8 + 6243748742141230639968*t**7 + 27761407182086143225024*t**6 + 70066148869420956398592*t**5 + 109701225644313784229376*t**4 + 109009005495588442152960*t**3 + 67072101084384786432000*t**2 + 23339979742629593088000*t + 3513592776846090240000
    ]

    assert groebner(I, x, t, order='grlex') == [
        16996618586000601590732959134095643086442*t**3*x - 32936701459297092865176560282688198064839*t**3 + 78592411049800639484139414821529525782364*t**2*x - 120753953358671750165454009478961405619916*t**2 + 120988399875140799712152158915653654637280*t*x - 144576390266626470824138354942076045758736*t + 60017634054270480831259316163620768960*x**2 + 61976058033571109604821862786675242894400*x - 56266268491293858791834120380427754600960,
        576689018321912327136790519059646508441672750656050290242749*t**4 + 2326673103677477425562248201573604572527893938459296513327336*t**3 + 110743790416688497407826310048520299245819959064297990236000*t**2*x + 3308669114229100853338245486174247752683277925010505284338016*t**2 + 323150205645687941261103426627818874426097912639158572428800*t*x + 1914335199925152083917206349978534224695445819017286960055680*t + 861662882561803377986838989464278045397192862768588480000*x**2 + 235296483281783440197069672204341465480107019878814196672000*x + 361850798943225141738895123621685122544503614946436727532800,
        -117584925286448670474763406733005510014188341867*t**3 + 68566565876066068463853874568722190223721653044*t**2*x - 435970731348366266878180788833437896139920683940*t**2 + 196297602447033751918195568051376792491869233408*t*x - 525011527660010557871349062870980202067479780112*t + 517905853447200553360289634770487684447317120*x**3 + 569119014870778921949288951688799397569321920*x**2 + 138877356748142786670127389526667463202210102080*x - 205109210539096046121625447192779783475018619520,
        -3725142681462373002731339445216700112264527*t**3 + 583711207282060457652784180668273817487940*t**2*x - 12381382393074485225164741437227437062814908*t**2 + 151081054097783125250959636747516827435040*t*x**2 + 1814103857455163948531448580501928933873280*t*x - 13353115629395094645843682074271212731433648*t + 236415091385250007660606958022544983766080*x**2 + 1390443278862804663728298060085399578417600*x - 4716885828494075789338754454248931750698880
    ]

@XFAIL
def test_benchmark_czichowski():
    skip('This takes too much time (without gmpy)')

    config.setup('GB_METHOD', 'f5b')
    helper_test_benchmark_czichowski()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_benchmark_czichowski()

def helper_test_benchmark_cyclic_4():
    a, b, c, d = symbols('a b c d')

    I = [a + b + c + d, a*b + a*d + b*c + b*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1]

    assert groebner(I, a, b, c, d, order='lex') == [
        4*a + 3*d**9 - 4*d**5 - 3*d,
        4*b + 4*c - 3*d**9 + 4*d**5 + 7*d,
        4*c**2 + 3*d**10 - 4*d**6 - 3*d**2,
        4*c*d**4 + 4*c - d**9 + 4*d**5 + 5*d, d**12 - d**8 - d**4 + 1
    ]

    assert groebner(I, a, b, c, d, order='grlex') == [
        3*b*c - c**2 + d**6 - 3*d**2,
        -b + 3*c**2*d**3 - c - d**5 - 4*d,
        -b + 3*c*d**4 + 2*c + 2*d**5 + 2*d,
        c**4 + 2*c**2*d**2 - d**4 - 2,
        c**3*d + c*d**3 + d**4 + 1,
        b*c**2 - c**3 - c**2*d - 2*c*d**2 - d**3,
        b**2 - c**2, b*d + c**2 + c*d + d**2,
        a + b + c + d
    ]

def test_benchmark_cyclic_4():
    config.setup('GB_METHOD', 'f5b')
    helper_test_benchmark_cyclic_4()
    config.setup('GB_METHOD', 'buchberger')
    helper_test_benchmark_cyclic_4()

def test_sig_key():
    s1 = sig((0,) * 3, 2)
    s2 = sig((1,) * 3, 4)
    s3 = sig((2,) * 3, 2)

    assert sig_key(s1, lex) > sig_key(s2, lex)
    assert sig_key(s2, lex) < sig_key(s3, lex)

def test_lbp_key():
    p1 = lbp(sig((0,) * 4, 3), [], 12)
    p2 = lbp(sig((0,) * 4, 4), [], 13)
    p3 = lbp(sig((0,) * 4, 4), [], 12)

    assert lbp_key(p1, lex) > lbp_key(p2, lex)
    assert lbp_key(p2, lex) < lbp_key(p3, lex)

def test_critical_pair():
    # from cyclic4 with grlex
    p1 = (((0, 0, 0, 0), 4), [((0, 1, 1, 2), QQ(1,1)), ((0, 0, 2, 2), QQ(1,1)), ((0, 0, 0, 4), QQ(-1,1)), ((0, 0, 0, 0), QQ(-1,1))], 4)
    q1 = (((0, 0, 0, 0), 2), [((0, 2, 0, 0), QQ(-1,1)), ((0, 1, 0, 1), QQ(-1,1)), ((0, 0, 1, 1), QQ(-1,1)), ((0, 0, 0, 2), QQ(-1,1))], 2)

    p2 = (((0, 0, 0, 2), 3), [((0, 0, 3, 2), QQ(1,1)), ((0, 0, 2, 3), QQ(1,1)), ((0, 0, 1, 0), QQ(-1,1)), ((0, 0, 0, 1), QQ(-1,1))], 5)
    q2 = (((0, 0, 2, 2), 2), [((0, 0, 1, 5), QQ(1,1)), ((0, 0, 0, 6), QQ(1,1)), ((0, 1, 1, 0), QQ(1,1)), ((0, 0, 1, 1), QQ(1,1))], 13)

    assert critical_pair(p1, q1, 3, grlex, QQ) == (((0, 0, 1, 2), 2), ((0, 0, 1, 2), QQ(-1,1)), (((0, 0, 0, 0), 2), [((0, 2, 0, 0), QQ(-1,1)), ((0, 1, 0, 1), QQ(-1,1)), ((0, 0, 1, 1), QQ(-1,1)), ((0, 0, 0, 2), QQ(-1,1))], 2), ((0, 1, 0, 0), 4), ((0, 1, 0, 0), QQ(1,1)), (((0, 0, 0, 0), 4), [((0, 1, 1, 2), QQ(1,1)), ((0, 0, 2, 2), QQ(1,1)), ((0, 0, 0, 4), QQ(-1,1)), ((0, 0, 0, 0), QQ(-1,1))], 4))
    assert critical_pair(p2, q2, 3, grlex, QQ) == (((0, 0, 4, 2), 2), ((0, 0, 2, 0), QQ(1,1)), (((0, 0, 2, 2), 2), [((0, 0, 1, 5), QQ(1,1)), ((0, 0, 0, 6), QQ(1,1)), ((0, 1, 1, 0), QQ(1,1)), ((0, 0, 1, 1), QQ(1,1))], 13), ((0, 0, 0, 5), 3), ((0, 0, 0, 3), QQ(1,1)), (((0, 0, 0, 2), 3), [((0, 0, 3, 2), QQ(1,1)), ((0, 0, 2, 3), QQ(1,1)), ((0, 0, 1, 0), QQ(-1,1)), ((0, 0, 0, 1), QQ(-1,1))], 5))

def test_cp_key():
    # from cyclic4 with grlex
    p1 = (((0, 0, 0, 0), 4), [((0, 1, 1, 2), QQ(1,1)), ((0, 0, 2, 2), QQ(1,1)), ((0, 0, 0, 4), QQ(-1,1)), ((0, 0, 0, 0), QQ(-1,1))], 4)
    q1 = (((0, 0, 0, 0), 2), [((0, 2, 0, 0), QQ(-1,1)), ((0, 1, 0, 1), QQ(-1,1)), ((0, 0, 1, 1), QQ(-1,1)), ((0, 0, 0, 2), QQ(-1,1))], 2)

    p2 = (((0, 0, 0, 2), 3), [((0, 0, 3, 2), QQ(1,1)), ((0, 0, 2, 3), QQ(1,1)), ((0, 0, 1, 0), QQ(-1,1)), ((0, 0, 0, 1), QQ(-1,1))], 5)
    q2 = (((0, 0, 2, 2), 2), [((0, 0, 1, 5), QQ(1,1)), ((0, 0, 0, 6), QQ(1,1)), ((0, 1, 1, 0), QQ(1,1)), ((0, 0, 1, 1), QQ(1,1))], 13)

    cp1 = critical_pair(p1, q1, 3, grlex, QQ)
    cp2 = critical_pair(p2, q2, 3, grlex, QQ)

    assert cp_key(cp1, grlex) < cp_key(cp2, grlex)

    cp1 = critical_pair(p1, p2, 3, grlex, QQ)
    cp2 = critical_pair(q1, q2, 3, grlex, QQ)

    assert cp_key(cp1, grlex) < cp_key(cp2, grlex)

def test_is_rewritable_or_comparable():
    # from katsura4 with grlex
    p = lbp(sig((0, 0, 2, 1), 2), [], 2)
    B = [lbp(sig((0, 0, 0, 1), 2), [((0, 0, 2, 1), QQ(1,1)), ((0, 0, 1, 2), QQ(76,35)), ((0, 0, 0, 3), QQ(13,7)), ((0, 2, 0, 0), QQ(2,45)), ((0, 1, 1, 0), QQ(1,5)), ((0, 1, 0, 1), QQ(5,63)), ((0, 0, 2, 0), QQ(4,45)), ((0, 0, 1, 1), QQ(-32,105)), ((0, 0, 0, 2), QQ(-13,21))], 6)]

    # rewritable:
    assert is_rewritable_or_comparable(Sign(p), Num(p), B, 3, QQ) == True

    p = lbp(sig((0, 1, 1, 0), 2), [], 7)
    B = [lbp(sig((0, 0, 0, 0), 3), [((0, 1, 1, 0), QQ(10,3)), ((0, 1, 0, 1), QQ(4,3)), ((0, 0, 2, 0), QQ(4,1)), ((0, 0, 1, 1), QQ(22,3)), ((0, 0, 0, 2), QQ(4,1)), ((0, 1, 0, 0), QQ(-1,3)), ((0, 0, 1, 0), QQ(-4,3)), ((0, 0, 0, 1), QQ(-4,3))], 3)]
    # comparable:
    assert is_rewritable_or_comparable(Sign(p), Num(p), B, 3, QQ) == True

def test_f5_reduce():
    # katsura3 with lex
    F = [(((0, 0, 0), 1), [((1, 0, 0), QQ(1,1)), ((0, 1, 0), QQ(2,1)), ((0, 0, 1), QQ(2,1)), ((0, 0, 0), QQ(-1,1))], 1), (((0, 0, 0), 2), [((0, 2, 0), QQ(6,1)), ((0, 1, 1), QQ(8,1)), ((0, 1, 0), QQ(-2,1)), ((0, 0, 2), QQ(6,1)), ((0, 0, 1), QQ(-2,1))], 2), (((0, 0, 0), 3), [((0, 1, 1), QQ(10,3)), ((0, 1, 0), QQ(-1,3)), ((0, 0, 2), QQ(4,1)), ((0, 0, 1), QQ(-4,3))], 3), (((0, 0, 1), 2), [((0, 1, 0), QQ(1,1)), ((0, 0, 3), QQ(30,1)), ((0, 0, 2), QQ(-79,7)), ((0, 0, 1), QQ(3,7))], 4), (((0, 0, 2), 2), [((0, 0, 4), QQ(1,1)), ((0, 0, 3), QQ(-10,21)), ((0, 0, 2), QQ(1,84)), ((0, 0, 1), QQ(1,84))], 5)]

    cp = critical_pair(F[0], F[1], 2, lex, QQ)
    s = s_poly(cp, 2, lex, QQ)

    assert f5_reduce(s, F, 2, lex, QQ) == (((0, 2, 0), 1), [], 1)

    s = lbp(sig(Sign(s)[0], 100), Polyn(s), Num(s))
    assert f5_reduce(s, F, 2, lex, QQ) == s

def test_matrix_fglm():
    pass  # see test_polytools.py

def test_representing_matrices():
    basis = [(0, 0), (0, 1), (1, 0), (1, 1)]
    F = [[((2, 0), QQ(1,1)), ((1, 0), QQ(-1,1)), ((0, 1), QQ(-3,1)), ((0, 0), QQ(1,1))],
        [((0, 2), QQ(1,1)), ((1, 0), QQ(-2,1)), ((0, 1), QQ(1,1)), ((0, 0), QQ(-1,1))]]

    assert _representing_matrices(basis, F, 1, grlex, QQ) ==[ \
        [[QQ(0,1), QQ(0,1), QQ(-1,1), QQ(3,1)],
        [QQ(0,1), QQ(0,1), QQ(3,1), QQ(-4,1)],
        [QQ(1,1), QQ(0,1), QQ(1,1), QQ(6,1)],
        [QQ(0,1), QQ(1,1), QQ(0,1), QQ(1,1)]],
        [[QQ(0,1), QQ(1,1), QQ(0,1), QQ(-2,1)],
        [QQ(1,1), QQ(-1,1), QQ(0,1), QQ(6,1)],
        [QQ(0,1), QQ(2,1), QQ(0,1), QQ(3,1)],
        [QQ(0,1), QQ(0,1), QQ(1,1), QQ(-1,1)]]]
