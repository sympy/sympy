from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.modulargcd import modgcd_univariate

def test_modgcd_univariate_integers():
    R, x = ring("x", ZZ)

    f, g = R.zero, R.zero
    assert modgcd_univariate(f, g) == (0, 0, 0)

    f, g = R.zero, x
    assert modgcd_univariate(f, g) == (x, 0, 1)
    assert modgcd_univariate(g, f) == (x, 1, 0)

    f, g = R.zero, -x
    assert modgcd_univariate(f, g) == (x, 0, -1)
    assert modgcd_univariate(g, f) == (x, -1, 0)

    f, g = 2*x, R(2)
    assert modgcd_univariate(f, g) == (2, x, 1)

    f, g = 2*x + 2, 6*x**2 - 6
    assert modgcd_univariate(f, g) == (2*x + 2, 1, 3*x - 3)

    f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
    g = x**3 + 6*x**2 + 11*x + 6

    h = x**2 + 3*x + 2

    cff = x**2 + 5*x + 4
    cfg = x + 3

    assert modgcd_univariate(f, g) == (h, cff, cfg)

    f = x**4 - 4
    g = x**4 + 4*x**2 + 4

    h = x**2 + 2

    cff = x**2 - 2
    cfg = x**2 + 2

    assert modgcd_univariate(f, g) == (h, cff, cfg)

    f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    h = 1

    cff = f
    cfg = g

    assert modgcd_univariate(f, g) == (h, cff, cfg)

    f = - 352518131239247345597970242177235495263669787845475025293906825864749649589178600387510272*x**49 \
        + 46818041807522713962450042363465092040687472354933295397472942006618953623327997952*x**42 \
        + 378182690892293941192071663536490788434899030680411695933646320291525827756032*x**35 \
        + 112806468807371824947796775491032386836656074179286744191026149539708928*x**28 \
        - 12278371209708240950316872681744825481125965781519138077173235712*x**21 \
        + 289127344604779611146960547954288113529690984687482920704*x**14 \
        + 19007977035740498977629742919480623972236450681*x**7 \
        + 311973482284542371301330321821976049

    g =   365431878023781158602430064717380211405897160759702125019136*x**21 \
        + 197599133478719444145775798221171663643171734081650688*x**14 \
        - 9504116979659010018253915765478924103928886144*x**7 \
        - 311973482284542371301330321821976049

    assert modgcd_univariate(f, f.diff(x))[0] == g

    f = 1317378933230047068160*x + 2945748836994210856960
    g = 120352542776360960*x + 269116466014453760

    h = 120352542776360960*x + 269116466014453760
    cff = 10946
    cfg = 1

    assert modgcd_univariate(f, g) == (h, cff, cfg)
