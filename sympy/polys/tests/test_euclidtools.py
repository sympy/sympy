"""Tests for Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences. """

from sympy.polys.euclidtools import (
    dup_gcdex, dup_half_gcdex, dup_invert,
    dup_euclidean_prs, dmp_euclidean_prs,
    dup_primitive_prs, dmp_primitive_prs,
    dup_subresultants, dmp_subresultants,
    dup_prs_resultant, dmp_prs_resultant,
    dmp_zz_collins_resultant,
    dmp_qq_collins_resultant,
    dup_resultant, dmp_resultant,
    dup_discriminant, dmp_discriminant,
    dup_zz_heu_gcd, dmp_zz_heu_gcd,
    dup_qq_heu_gcd, dmp_qq_heu_gcd,
    dup_rr_prs_gcd, dmp_rr_prs_gcd,
    dup_ff_prs_gcd, dmp_ff_prs_gcd,
    dup_inner_gcd, dmp_inner_gcd,
    dup_lcm, dmp_lcm,
    dmp_content, dmp_primitive)

from sympy.polys.densebasic import (
    dmp_one_p,
    dup_LC, dmp_LC,
    dup_normal, dmp_normal)

from sympy.polys.densearith import (
    dup_add,
    dup_mul, dmp_mul,
    dup_exquo)

from sympy.polys.densetools import (
    dup_diff)

from sympy.polys.specialpolys import (
    f_4, f_5, f_6,
    dmp_fateman_poly_F_1,
    dmp_fateman_poly_F_2,
    dmp_fateman_poly_F_3)

from sympy.polys.domains import ZZ, QQ

def test_dup_gcdex():
    f = dup_normal([1,-2,-6,12,15], QQ)
    g = dup_normal([1,1,-4,-4], QQ)

    s = [QQ(-1,5),QQ(3,5)]
    t = [QQ(1,5),QQ(-6,5),QQ(2)]
    h = [QQ(1),QQ(1)]

    assert dup_half_gcdex(f, g, QQ) == (s, h)
    assert dup_gcdex(f, g, QQ) == (s, t, h)

    f = dup_normal([1,4,0,-1,1], QQ)
    g = dup_normal([1,0,-1,1], QQ)

    s, t, h = dup_gcdex(f, g, QQ)
    S, T, H = dup_gcdex(g, f, QQ)

    assert dup_add(dup_mul(s, f, QQ),
                   dup_mul(t, g, QQ), QQ) == h
    assert dup_add(dup_mul(S, g, QQ),
                   dup_mul(T, f, QQ), QQ) == H

    f = dup_normal([2,0], QQ)
    g = dup_normal([1,0,-16], QQ)

    s = [QQ(1,32),QQ(0)]
    t = [QQ(-1,16)]
    h = [QQ(1)]

    assert dup_half_gcdex(f, g, QQ) == (s, h)
    assert dup_gcdex(f, g, QQ) == (s, t, h)

def test_dup_invert():
    assert dup_invert([QQ(2),QQ(0)], [QQ(1),QQ(0),QQ(-16)], QQ) == [QQ(1,32),QQ(0)]

def test_dup_euclidean_prs():
    f = QQ.map([1, 0, 1, 0, -3, -3, 8, 2, -5])
    g = QQ.map([3, 0, 5, 0, -4, -9, 21])

    assert dup_euclidean_prs(f, g, QQ) == [f, g,
        [-QQ(5,9), QQ(0,1), QQ(1,9), QQ(0,1), -QQ(1,3)],
        [-QQ(117,25), -QQ(9,1), QQ(441,25)],
        [QQ(233150,19773), -QQ(102500,6591)],
        [-QQ(1288744821,543589225)]]

def test_dup_primitive_prs():
    f = ZZ.map([1, 0, 1, 0, -3, -3, 8, 2, -5])
    g = ZZ.map([3, 0, 5, 0, -4, -9, 21])

    assert dup_primitive_prs(f, g, ZZ) == [f, g,
        [-ZZ(5), ZZ(0), ZZ(1), ZZ(0), -ZZ(3)],
        [ZZ(13), ZZ(25), -ZZ(49)],
        [ZZ(4663), -ZZ(6150)],
        [ZZ(1)]]

def test_dup_subresultants():
    assert dup_resultant([], [], ZZ) == ZZ(0)

    assert dup_resultant([ZZ(1)], [], ZZ) == ZZ(0)
    assert dup_resultant([], [ZZ(1)], ZZ) == ZZ(0)

    f = dup_normal([1,0,1,0,-3,-3,8,2,-5], ZZ)
    g = dup_normal([3,0,5,0,-4,-9,21], ZZ)

    a = dup_normal([15,0,-3,0,9], ZZ)
    b = dup_normal([65,125,-245], ZZ)
    c = dup_normal([9326,-12300], ZZ)
    d = dup_normal([260708], ZZ)

    assert dup_subresultants(f, g, ZZ) == [f, g, a, b, c, d]
    assert dup_resultant(f, g, ZZ) == dup_LC(d, ZZ)

    f = dup_normal([1,-2,1], ZZ)
    g = dup_normal([1,0,-1], ZZ)

    a = dup_normal([2,-2], ZZ)

    assert dup_subresultants(f, g, ZZ) == [f, g, a]
    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([1,0, 1], ZZ)
    g = dup_normal([1,0,-1], ZZ)

    a = dup_normal([-2], ZZ)

    assert dup_subresultants(f, g, ZZ) ==  [f, g, a]
    assert dup_resultant(f, g, ZZ) == 4

    f = dup_normal([1,0,-1], ZZ)
    g = dup_normal([1,-1,0,2], ZZ)

    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([3,0,-1,0], ZZ)
    g = dup_normal([5,0,1], ZZ)

    assert dup_resultant(f, g, ZZ) == 64

    f = dup_normal([1,-2,7], ZZ)
    g = dup_normal([1,0,-1,5], ZZ)

    assert dup_resultant(f, g, ZZ) == 265

    f = dup_normal([1,-6,11,-6], ZZ)
    g = dup_normal([1,-15,74,-120], ZZ)

    assert dup_resultant(f, g, ZZ) == -8640

    f = dup_normal([1,-6,11,-6], ZZ)
    g = dup_normal([1,-10,29,-20], ZZ)

    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([1,0,0,-1], ZZ)
    g = dup_normal([1,2,2,-1], ZZ)

    assert dup_resultant(f, g, ZZ) == 16

    f = dup_normal([1,0,0,0,0,0,0,0,-2], ZZ)
    g = dup_normal([1,-1], ZZ)

    assert dup_resultant(f, g, ZZ) == -1

def test_dmp_subresultants():
    assert dmp_resultant([[]], [[]], 1, ZZ) == []
    assert dmp_prs_resultant([[]], [[]], 1, ZZ)[0] == []
    assert dmp_zz_collins_resultant([[]], [[]], 1, ZZ) == []
    assert dmp_qq_collins_resultant([[]], [[]], 1, ZZ) == []

    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []
    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []
    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []

    assert dmp_resultant([[]], [[ZZ(1)]], 1, ZZ) == []
    assert dmp_prs_resultant([[]], [[ZZ(1)]], 1, ZZ)[0] == []
    assert dmp_zz_collins_resultant([[]], [[ZZ(1)]], 1, ZZ) == []
    assert dmp_qq_collins_resultant([[]], [[ZZ(1)]], 1, ZZ) == []

    f = dmp_normal([[3,0],[],[-1,0,0,-4]], 1, ZZ)
    g = dmp_normal([[1],[1,0,0,0],[-9]], 1, ZZ)

    a = dmp_normal([[3,0,0,0,0],[1,0,-27,4]], 1, ZZ)
    b = dmp_normal([[-3,0,0,-12,1,0,-54,8,729,-216,16]], 1, ZZ)

    r = dmp_LC(b, ZZ)

    assert dmp_subresultants(f, g, 1, ZZ) == [f, g, a, b]

    assert dmp_resultant(f, g, 1, ZZ) == r
    assert dmp_prs_resultant(f, g, 1, ZZ)[0] == r
    assert dmp_zz_collins_resultant(f, g, 1, ZZ) == r
    assert dmp_qq_collins_resultant(f, g, 1, ZZ) == r

    f = dmp_normal([[-1],[],[],[5]], 1, ZZ)
    g = dmp_normal([[3,1],[],[]], 1, ZZ)

    a = dmp_normal([[45,30,5]], 1, ZZ)
    b = dmp_normal([[675,675,225,25]], 1, ZZ)

    r = dmp_LC(b, ZZ)

    assert dmp_subresultants(f, g, 1, ZZ) == [f, g, a]
    assert dmp_resultant(f, g, 1, ZZ) == r
    assert dmp_prs_resultant(f, g, 1, ZZ)[0] == r
    assert dmp_zz_collins_resultant(f, g, 1, ZZ) == r
    assert dmp_qq_collins_resultant(f, g, 1, ZZ) == r

    f = [[[[[6]]]], [[[[-3]]], [[[-2]], [[]]]], [[[[1]], [[]]], [[[]]]]]
    g = [[[[[1]]]], [[[[-1], [-1, 0]]]], [[[[1, 0], []]]]]

    r = [[[[1]], [[-3], [-3, 0]], [[9, 0], []]], [[[-2], [-2, 0]], [[6],
         [12, 0], [6, 0, 0]], [[-18, 0], [-18, 0, 0], []]], [[[4, 0],
         []], [[-12, 0], [-12, 0, 0], []], [[36, 0, 0], [], []]]]

    assert dmp_zz_collins_resultant(f, g, 4, ZZ) == r

    f = [[[[[QQ(1,1)]]]], [[[[QQ(-1,2)]]], [[[QQ(-1,3)]], [[]]]], [[[[QQ(1,6)]], [[]]], [[[]]]]]
    g = [[[[[QQ(1,1)]]]], [[[[QQ(-1,1)], [QQ(-1,1), QQ(0, 1)]]]], [[[[QQ(1,1), QQ(0,1)], []]]]]

    r = [[[[QQ(1,36)]], [[QQ(-1,12)], [QQ(-1,12), QQ(0,1)]], [[QQ(1,4), QQ(0,1)], []]],
         [[[QQ(-1,18)], [QQ(-1,18), QQ(0,1)]], [[QQ(1,6)], [QQ(1,3), QQ(0,1)], [QQ(1,6),
            QQ(0,1), QQ(0,1)]], [[QQ(-1,2), QQ(0,1)], [QQ(-1,2), QQ(0,1), QQ(0,1)], []]],
         [[[QQ(1,9), QQ(0,1)], []], [[QQ(-1,3), QQ(0,1)], [QQ(-1,3), QQ(0,1), QQ(0,1)], []],
          [[QQ(1,1), QQ(0,1), QQ(0,1)], [], []]]]

    assert dmp_qq_collins_resultant(f, g, 4, QQ) == r

def test_dup_discriminant():
    assert dup_discriminant([], ZZ) == 0
    assert dup_discriminant([1,0], ZZ) == 1

    assert dup_discriminant([1,3,9,-13], ZZ) == -11664
    assert dup_discriminant([5,0,1,0,0,2], ZZ) == 31252160
    assert dup_discriminant([1,2,6,-22,13], ZZ) == 0
    assert dup_discriminant([12,0,0,15,30,1,0,1], ZZ) == -220289699947514112

def test_dmp_discriminant():
    assert dmp_discriminant([], 0, ZZ) == 0
    assert dmp_discriminant([[]], 1, ZZ) == []

    assert dmp_discriminant([[1,0]], 1, ZZ) == []

    assert dmp_discriminant([1,3,9,-13], 0, ZZ) == -11664
    assert dmp_discriminant([5,0,1,0,0,2], 0, ZZ) == 31252160
    assert dmp_discriminant([1,2,6,-22,13], 0, ZZ) == 0
    assert dmp_discriminant([12,0,0,15,30,1,0,1], 0, ZZ) == -220289699947514112

    assert dmp_discriminant([[1,0],[],[2,0]], 1, ZZ) == [-8,0,0]
    assert dmp_discriminant([[1,0,2],[]], 1, ZZ) == [1]

    assert dmp_discriminant([[[1],[]],[[1,0]]], 2, ZZ) == [[1]]

    assert dmp_discriminant([[[[1]],[[]]],[[[1],[]]],[[[1,0]]]], 3, ZZ) == \
        [[[-4, 0]], [[1], [], []]]
    assert dmp_discriminant([[[[[1]]],[[[]]]],[[[[1]],[[]]]],[[[[1],[]]]],[[[[1,0]]]]], 4, ZZ) == \
        [[[[-27,0,0]]],[[[18,0],[]],[[-4],[],[],[]]],[[[-4,0]],[[1],[],[]],[[]],[[]]]]

def test_dup_gcd():
    assert dup_zz_heu_gcd([], [], ZZ) == ([], [], [])
    assert dup_rr_prs_gcd([], [], ZZ) == ([], [], [])

    assert dup_zz_heu_gcd([2], [], ZZ) == ([2], [1], [])
    assert dup_rr_prs_gcd([2], [], ZZ) == ([2], [1], [])

    assert dup_zz_heu_gcd([-2], [], ZZ) == ([2], [-1], [])
    assert dup_rr_prs_gcd([-2], [], ZZ) == ([2], [-1], [])

    assert dup_zz_heu_gcd([], [-2], ZZ) == ([2], [], [-1])
    assert dup_rr_prs_gcd([], [-2], ZZ) == ([2], [], [-1])

    assert dup_zz_heu_gcd([], [2,4], ZZ) == ([2,4], [], [1])
    assert dup_rr_prs_gcd([], [2,4], ZZ) == ([2,4], [], [1])

    assert dup_zz_heu_gcd([2,4], [], ZZ) == ([2,4], [1], [])
    assert dup_rr_prs_gcd([2,4], [], ZZ) == ([2,4], [1], [])

    assert dup_zz_heu_gcd([2], [2], ZZ) == ([2], [1], [1])
    assert dup_rr_prs_gcd([2], [2], ZZ) == ([2], [1], [1])

    assert dup_zz_heu_gcd([-2], [2], ZZ) == ([2], [-1], [1])
    assert dup_rr_prs_gcd([-2], [2], ZZ) == ([2], [-1], [1])

    assert dup_zz_heu_gcd([2], [-2], ZZ) == ([2], [1], [-1])
    assert dup_rr_prs_gcd([2], [-2], ZZ) == ([2], [1], [-1])

    assert dup_zz_heu_gcd([-2], [-2], ZZ) == ([2], [-1], [-1])
    assert dup_rr_prs_gcd([-2], [-2], ZZ) == ([2], [-1], [-1])

    assert dup_zz_heu_gcd([1,2,1], [1], ZZ) == ([1], [1, 2, 1], [1])
    assert dup_rr_prs_gcd([1,2,1], [1], ZZ) == ([1], [1, 2, 1], [1])

    assert dup_zz_heu_gcd([1,2,1], [2], ZZ) == ([1], [1, 2, 1], [2])
    assert dup_rr_prs_gcd([1,2,1], [2], ZZ) == ([1], [1, 2, 1], [2])

    assert dup_zz_heu_gcd([2,4,2], [2], ZZ) == ([2], [1, 2, 1], [1])
    assert dup_rr_prs_gcd([2,4,2], [2], ZZ) == ([2], [1, 2, 1], [1])

    assert dup_zz_heu_gcd([2], [2,4,2], ZZ) == ([2], [1], [1, 2, 1])
    assert dup_rr_prs_gcd([2], [2,4,2], ZZ) == ([2], [1], [1, 2, 1])

    assert dup_zz_heu_gcd([2,4,2], [1,1], ZZ) == ([1, 1], [2, 2], [1])
    assert dup_rr_prs_gcd([2,4,2], [1,1], ZZ) == ([1, 1], [2, 2], [1])

    assert dup_zz_heu_gcd([1,1], [2,4,2], ZZ) == ([1, 1], [1], [2, 2])
    assert dup_rr_prs_gcd([1,1], [2,4,2], ZZ) == ([1, 1], [1], [2, 2])

    f, g = [1, -31], [1, 0]

    assert dup_zz_heu_gcd(f, g, ZZ) == ([1], f, g)
    assert dup_rr_prs_gcd(f, g, ZZ) == ([1], f, g)

    f = [1,8,21,22,8]
    g = [1,6,11,6]

    h = [1,3,2]

    cff = [1,5,4]
    cfg = [1,3]

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = [1,0,0,0,-4]
    g = [1,0,4,0, 4]

    h = [1,0,2]

    cff = [1,0,-2]
    cfg = [1,0, 2]

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = [1,0,1,0,-3,-3,8,2,-5]
    g = [3,0,5,-0,-4,-9,21]

    h = [1]

    cff = f
    cfg = g

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = dup_normal([1,0,1,0,-3,-3,8,2,-5], QQ)
    g = dup_normal([3,0,5,-0,-4,-9,21], QQ)

    h = dup_normal([1], QQ)

    assert dup_qq_heu_gcd(f, g, QQ) == (h, cff, cfg)
    assert dup_ff_prs_gcd(f, g, QQ) == (h, cff, cfg)

    f = [-352518131239247345597970242177235495263669787845475025293906825864749649589178600387510272,
         0, 0, 0, 0, 0, 0,
         46818041807522713962450042363465092040687472354933295397472942006618953623327997952,
         0, 0, 0, 0, 0, 0,
         378182690892293941192071663536490788434899030680411695933646320291525827756032,
         0, 0, 0, 0, 0, 0,
         112806468807371824947796775491032386836656074179286744191026149539708928,
         0, 0, 0, 0, 0, 0,
         -12278371209708240950316872681744825481125965781519138077173235712,
         0, 0, 0, 0, 0, 0,
         289127344604779611146960547954288113529690984687482920704,
         0, 0, 0, 0, 0, 0,
         19007977035740498977629742919480623972236450681,
         0, 0, 0, 0, 0, 0,
         311973482284542371301330321821976049]

    g = [365431878023781158602430064717380211405897160759702125019136,
         0, 0, 0, 0, 0, 0,
         197599133478719444145775798221171663643171734081650688,
         0, 0, 0, 0, 0, 0,
         -9504116979659010018253915765478924103928886144,
         0, 0, 0, 0, 0, 0,
         -311973482284542371301330321821976049]

    f = dup_normal(f, ZZ)
    g = dup_normal(g, ZZ)

    assert dup_zz_heu_gcd(f, dup_diff(f, 1, ZZ), ZZ)[0] == g
    assert dup_rr_prs_gcd(f, dup_diff(f, 1, ZZ), ZZ)[0] == g

    f = [QQ(1,2),QQ(1),QQ(1,2)]
    g = [QQ(1,2),QQ(1,2)]

    h = [QQ(1), QQ(1)]

    assert dup_qq_heu_gcd(f, g, QQ) == (h, g, [QQ(1,2)])
    assert dup_ff_prs_gcd(f, g, QQ) == (h, g, [QQ(1,2)])

def test_dmp_gcd():
    assert dmp_zz_heu_gcd([[]], [[]], 1, ZZ) == ([[]], [[]], [[]])
    assert dmp_rr_prs_gcd([[]], [[]], 1, ZZ) == ([[]], [[]], [[]])

    assert dmp_zz_heu_gcd([[2]], [[]], 1, ZZ) == ([[2]], [[1]], [[]])
    assert dmp_rr_prs_gcd([[2]], [[]], 1, ZZ) == ([[2]], [[1]], [[]])

    assert dmp_zz_heu_gcd([[-2]], [[]], 1, ZZ) == ([[2]], [[-1]], [[]])
    assert dmp_rr_prs_gcd([[-2]], [[]], 1, ZZ) == ([[2]], [[-1]], [[]])

    assert dmp_zz_heu_gcd([[]], [[-2]], 1, ZZ) == ([[2]], [[]], [[-1]])
    assert dmp_rr_prs_gcd([[]], [[-2]], 1, ZZ) == ([[2]], [[]], [[-1]])

    assert dmp_zz_heu_gcd([[]], [[2],[4]], 1, ZZ) == ([[2],[4]], [[]], [[1]])
    assert dmp_rr_prs_gcd([[]], [[2],[4]], 1, ZZ) == ([[2],[4]], [[]], [[1]])

    assert dmp_zz_heu_gcd([[2],[4]], [[]], 1, ZZ) == ([[2],[4]], [[1]], [[]])
    assert dmp_rr_prs_gcd([[2],[4]], [[]], 1, ZZ) == ([[2],[4]], [[1]], [[]])

    assert dmp_zz_heu_gcd([[2]], [[2]], 1, ZZ) == ([[2]], [[1]], [[1]])
    assert dmp_rr_prs_gcd([[2]], [[2]], 1, ZZ) == ([[2]], [[1]], [[1]])

    assert dmp_zz_heu_gcd([[-2]], [[2]], 1, ZZ) == ([[2]], [[-1]], [[1]])
    assert dmp_rr_prs_gcd([[-2]], [[2]], 1, ZZ) == ([[2]], [[-1]], [[1]])

    assert dmp_zz_heu_gcd([[2]], [[-2]], 1, ZZ) == ([[2]], [[1]], [[-1]])
    assert dmp_rr_prs_gcd([[2]], [[-2]], 1, ZZ) == ([[2]], [[1]], [[-1]])

    assert dmp_zz_heu_gcd([[-2]], [[-2]], 1, ZZ) == ([[2]], [[-1]], [[-1]])
    assert dmp_rr_prs_gcd([[-2]], [[-2]], 1, ZZ) == ([[2]], [[-1]], [[-1]])

    assert dmp_zz_heu_gcd([[1],[2],[1]], [[1]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[1]])
    assert dmp_rr_prs_gcd([[1],[2],[1]], [[1]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[1]])

    assert dmp_zz_heu_gcd([[1],[2],[1]], [[2]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[2]])
    assert dmp_rr_prs_gcd([[1],[2],[1]], [[2]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[2]])

    assert dmp_zz_heu_gcd([[2],[4],[2]], [[2]], 1, ZZ) == ([[2]], [[1], [2], [1]], [[1]])
    assert dmp_rr_prs_gcd([[2],[4],[2]], [[2]], 1, ZZ) == ([[2]], [[1], [2], [1]], [[1]])

    assert dmp_zz_heu_gcd([[2]], [[2],[4],[2]], 1, ZZ) == ([[2]], [[1]], [[1], [2], [1]])
    assert dmp_rr_prs_gcd([[2]], [[2],[4],[2]], 1, ZZ) == ([[2]], [[1]], [[1], [2], [1]])

    assert dmp_zz_heu_gcd([[2],[4],[2]], [[1],[1]], 1, ZZ) == ([[1], [1]], [[2], [2]], [[1]])
    assert dmp_rr_prs_gcd([[2],[4],[2]], [[1],[1]], 1, ZZ) == ([[1], [1]], [[2], [2]], [[1]])

    assert dmp_zz_heu_gcd([[1],[1]], [[2],[4],[2]], 1, ZZ) == ([[1], [1]], [[1]], [[2], [2]])
    assert dmp_rr_prs_gcd([[1],[1]], [[2],[4],[2]], 1, ZZ) == ([[1], [1]], [[1]], [[2], [2]])

    assert dmp_zz_heu_gcd([[[[1,2,1]]]], [[[[2,2]]]], 3, ZZ) == ([[[[1,1]]]], [[[[1,1]]]], [[[[2]]]])
    assert dmp_rr_prs_gcd([[[[1,2,1]]]], [[[[2,2]]]], 3, ZZ) == ([[[[1,1]]]], [[[[1,1]]]], [[[[2]]]])

    f, g = [[[[1,2,1],[1,1],[]]]], [[[[1,2,1]]]]
    h, cff, cfg = [[[[1,1]]]], [[[[1,1],[1],[]]]], [[[[1,1]]]]

    assert dmp_zz_heu_gcd(f, g, 3, ZZ) == (h, cff, cfg)
    assert dmp_rr_prs_gcd(f, g, 3, ZZ) == (h, cff, cfg)

    assert dmp_zz_heu_gcd(g, f, 3, ZZ) == (h, cfg, cff)
    assert dmp_rr_prs_gcd(g, f, 3, ZZ) == (h, cfg, cff)

    f, g, h = dmp_fateman_poly_F_1(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(4, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 4, ZZ)

    assert H == h and dmp_mul(H, cff, 4, ZZ) == f \
                  and dmp_mul(H, cfg, 4, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(6, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 6, ZZ)

    assert H == h and dmp_mul(H, cff, 6, ZZ) == f \
                  and dmp_mul(H, cfg, 6, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(8, ZZ)

    H, cff, cfg = dmp_zz_heu_gcd(f, g, 8, ZZ)

    assert H == h and dmp_mul(H, cff, 8, ZZ) == f \
                  and dmp_mul(H, cfg, 8, ZZ) == g

    f, g, h = dmp_fateman_poly_F_2(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_3(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_3(4, ZZ)
    H, cff, cfg = dmp_inner_gcd(f, g, 4, ZZ)

    assert H == h and dmp_mul(H, cff, 4, ZZ) == f \
                  and dmp_mul(H, cfg, 4, ZZ) == g

    f = [[QQ(1,2)],[QQ(1)],[QQ(1,2)]]
    g = [[QQ(1,2)],[QQ(1,2)]]

    h = [[QQ(1)],[QQ(1)]]

    assert dmp_qq_heu_gcd(f, g, 1, QQ) == (h, g, [[QQ(1,2)]])
    assert dmp_ff_prs_gcd(f, g, 1, QQ) == (h, g, [[QQ(1,2)]])

def test_dup_lcm():
    assert dup_lcm([2], [6], ZZ) == [6]

    assert dup_lcm([2,0,0,0], [6,0], ZZ) == [6,0,0,0]
    assert dup_lcm([2,0,0,0], [3,0], ZZ) == [6,0,0,0]

    assert dup_lcm([1,1,0], [1,0], ZZ) == [1,1,0]
    assert dup_lcm([1,1,0], [2,0], ZZ) == [2,2,0]
    assert dup_lcm([1,2,0], [1,0], ZZ) == [1,2,0]
    assert dup_lcm([2,1,0], [1,0], ZZ) == [2,1,0]
    assert dup_lcm([2,1,0], [2,0], ZZ) == [4,2,0]

def test_dmp_lcm():
    assert dmp_lcm([[2]], [[6]], 1, ZZ) == [[6]]
    assert dmp_lcm([[1],[]], [[1,0]], 1, ZZ) == [[1,0],[]]

    assert dmp_lcm([[2],[],[],[]], [[6,0,0],[]], 1, ZZ) == [[6,0,0],[],[],[]]
    assert dmp_lcm([[2],[],[],[]], [[3,0,0],[]], 1, ZZ) == [[6,0,0],[],[],[]]

    assert dmp_lcm([[1,0],[],[]], [[1,0,0],[]], 1, ZZ) == [[1,0,0],[],[]]

    f = [[2,-3,-2,3,0,0],[]]
    g = [[1,0,-2,0,1,0]]
    h = [[2,-3,-4,6,2,-3,0,0],[]]

    assert dmp_lcm(f, g, 1, ZZ) == h

    f = [[1],[-3,0],[-9,0,0],[-5,0,0,0]]
    g = [[1],[6,0],[12,0,0],[10,0,0,0],[3,0,0,0,0]]
    h = [[1],[1,0],[-18,0,0],[-50,0,0,0],[-47,0,0,0,0],[-15,0,0,0,0,0]]

    assert dmp_lcm(f, g, 1, ZZ) == h

def test_dmp_content():
    assert dmp_content([[-2]], 1, ZZ) == [2]

    f, g, F = [ZZ(3),ZZ(2),ZZ(1)], [ZZ(1)], []

    for i in xrange(0, 5):
        g = dup_mul(g, f, ZZ)
        F.insert(0, g)

    assert dmp_content(F, 1, ZZ) == f

    assert dmp_one_p(dmp_content(f_4, 2, ZZ), 1, ZZ)
    assert dmp_one_p(dmp_content(f_5, 2, ZZ), 1, ZZ)
    assert dmp_one_p(dmp_content(f_6, 3, ZZ), 2, ZZ)

def test_dmp_primitive():
    assert dmp_primitive([[]], 1, ZZ) == ([], [[]])
    assert dmp_primitive([[1]], 1, ZZ) == ([1], [[1]])

    f, g, F = [ZZ(3),ZZ(2),ZZ(1)], [ZZ(1)], []

    for i in xrange(0, 5):
        g = dup_mul(g, f, ZZ)
        F.insert(0, g)

    assert dmp_primitive(F, 1, ZZ) == (f,
        [ dup_exquo(c, f, ZZ) for c in F ])

    cont, f = dmp_primitive(f_4, 2, ZZ)
    assert dmp_one_p(cont, 1, ZZ) and f == f_4
    cont, f = dmp_primitive(f_5, 2, ZZ)
    assert dmp_one_p(cont, 1, ZZ) and f == f_5
    cont, f = dmp_primitive(f_6, 3, ZZ)
    assert dmp_one_p(cont, 2, ZZ) and f == f_6
