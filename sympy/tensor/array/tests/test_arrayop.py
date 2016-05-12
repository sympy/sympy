from sympy.utilities.pytest import raises

from sympy import symbols, sin, exp, log, cos, transpose, adjoint, conjugate
from sympy.tensor.array import Array

from sympy.tensor.array.arrayop import tensorproduct, tensorcontraction, derive_by_array, permutedims


def test_tensorproduct():
    x,y,z,t = symbols('x y z t')
    from sympy.abc import a,b,c,d
    assert tensorproduct() == 1
    assert tensorproduct([x]) == Array([x])
    assert tensorproduct([x], [y]) == Array([[x*y]])
    assert tensorproduct([x], [y], [z]) == Array([[[x*y*z]]])
    assert tensorproduct([x], [y], [z], [t]) == Array([[[[x*y*z*t]]]])

    assert tensorproduct(x) == x
    assert tensorproduct(x, y) == x*y
    assert tensorproduct(x, y, z) == x*y*z
    assert tensorproduct(x, y, z, t) == x*y*z*t

    A = Array([x, y])
    B = Array([1, 2, 3])
    C = Array([a, b, c, d])

    assert tensorproduct(A, B, C) == Array([[[a*x, b*x, c*x, d*x], [2*a*x, 2*b*x, 2*c*x, 2*d*x], [3*a*x, 3*b*x, 3*c*x, 3*d*x]],
                                            [[a*y, b*y, c*y, d*y], [2*a*y, 2*b*y, 2*c*y, 2*d*y], [3*a*y, 3*b*y, 3*c*y, 3*d*y]]])

    assert tensorproduct([x, y], [1, 2, 3]) == tensorproduct(A, B)

    assert tensorproduct(A, 2) == Array([2*x, 2*y])
    assert tensorproduct(A, [2]) == Array([[2*x], [2*y]])
    assert tensorproduct([2], A) == Array([[2*x, 2*y]])
    assert tensorproduct(a, A) == Array([a*x, a*y])
    assert tensorproduct(a, A, B) == Array([[a*x, 2*a*x, 3*a*x], [a*y, 2*a*y, 3*a*y]])
    assert tensorproduct(A, B, a) == Array([[a*x, 2*a*x, 3*a*x], [a*y, 2*a*y, 3*a*y]])
    assert tensorproduct(B, a, A) == Array([[a*x, a*y], [2*a*x, 2*a*y], [3*a*x, 3*a*y]])


def test_tensorcontraction():
    from sympy.abc import a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x
    B = Array(range(18), (2, 3, 3))
    assert tensorcontraction(B, (1, 2)) == Array([12, 39])
    C1 = Array([a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x], (2, 3, 2, 2))

    assert tensorcontraction(C1, (0, 2)) == Array([[a + o, b + p], [e + s, f + t], [i + w, j + x]])
    assert tensorcontraction(C1, (0, 2, 3)) == Array([a + p, e + t, i + x])
    assert tensorcontraction(C1, (2, 3)) == Array([[a + d, e + h, i + l], [m + p, q + t, u + x]])


def test_derivative_by_array():
    from sympy.abc import a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z

    bexpr = x*y**2*exp(z)*log(t)
    sexpr = sin(bexpr)
    cexpr = cos(bexpr)

    a = Array([sexpr])

    assert derive_by_array(sexpr, t) == x*y**2*exp(z)*cos(x*y**2*exp(z)*log(t))/t
    assert derive_by_array(sexpr, [x, y, z]) == Array([bexpr/x*cexpr, 2*y*bexpr/y**2*cexpr, bexpr*cexpr])
    assert derive_by_array(a, [x, y, z]) == Array([[bexpr/x*cexpr], [2*y*bexpr/y**2*cexpr], [bexpr*cexpr]])

    assert derive_by_array(sexpr, [[x, y], [z, t]]) == Array([[bexpr/x*cexpr, 2*y*bexpr/y**2*cexpr], [bexpr*cexpr, bexpr/log(t)/t*cexpr]])
    assert derive_by_array(a, [[x, y], [z, t]]) == Array([[[bexpr/x*cexpr], [2*y*bexpr/y**2*cexpr]], [[bexpr*cexpr], [bexpr/log(t)/t*cexpr]]])
    assert derive_by_array([[x, y], [z, t]], [x, y]) == Array([[[1, 0], [0, 0]], [[0, 1], [0, 0]]])
    assert derive_by_array([[x, y], [z, t]], [[x, y], [z, t]]) == Array([[[[1, 0], [0, 0]], [[0, 1], [0, 0]]],
                                                                         [[[0, 0], [1, 0]], [[0, 0], [0, 1]]]])


def test_issue_emerged_while_discussing_10972():
    ua = Array([-1,0])
    Fa = Array([[0, 1], [-1, 0]])
    po = tensorproduct(Fa, ua, Fa, ua)
    assert tensorcontraction(po, (1, 2), (4, 5)) == Array([[0, 0], [0, 1]])

    sa = symbols('a0:144')
    po = Array(sa, [2, 2, 3, 3, 2, 2])
    assert tensorcontraction(po, (0, 1), (2, 3), (4, 5)) == sa[0] + sa[108] + sa[111] + sa[124] + sa[127] + sa[140] + sa[143] + sa[16] + sa[19] + sa[3] + sa[32] + sa[35]
    assert tensorcontraction(po, (0, 1, 4, 5), (2, 3)) == sa[0] + sa[111] + sa[127] + sa[143] + sa[16] + sa[32]
    assert tensorcontraction(po, (0, 1), (4, 5)) == Array([[sa[0] + sa[108] + sa[111] + sa[3], sa[112] + sa[115] + sa[4] + sa[7],
                                                             sa[11] + sa[116] + sa[119] + sa[8]], [sa[12] + sa[120] + sa[123] + sa[15],
                                                             sa[124] + sa[127] + sa[16] + sa[19], sa[128] + sa[131] + sa[20] + sa[23]],
                                                            [sa[132] + sa[135] + sa[24] + sa[27], sa[136] + sa[139] + sa[28] + sa[31],
                                                             sa[140] + sa[143] + sa[32] + sa[35]]])
    assert tensorcontraction(po, (0, 1), (2, 3)) == Array([[sa[0] + sa[108] + sa[124] + sa[140] + sa[16] + sa[32], sa[1] + sa[109] + sa[125] + sa[141] + sa[17] + sa[33]],
                                                           [sa[110] + sa[126] + sa[142] + sa[18] + sa[2] + sa[34], sa[111] + sa[127] + sa[143] + sa[19] + sa[3] + sa[35]]])


def test_array_permutedims():
    sa = symbols('a0:144')

    m1 = Array(sa[:6], (2, 3))
    assert permutedims(m1, (1, 0)) == transpose(m1)
    assert m1.tomatrix().T == permutedims(m1, (1, 0)).tomatrix()

    assert m1.tomatrix().T == transpose(m1).tomatrix()
    assert m1.tomatrix().C == conjugate(m1).tomatrix()
    assert m1.tomatrix().H == adjoint(m1).tomatrix()

    assert m1.tomatrix().T == m1.transpose().tomatrix()
    assert m1.tomatrix().C == m1.conjugate().tomatrix()
    assert m1.tomatrix().H == m1.adjoint().tomatrix()

    po = Array(sa, [2, 2, 3, 3, 2, 2])

    raises(ValueError, lambda: permutedims(po, (1, 1)))
    raises(ValueError, lambda: po.transpose())
    raises(ValueError, lambda: po.adjoint())

    assert permutedims(po, reversed(range(po.rank()))) == Array(
        [[[[[[sa[0], sa[72]], [sa[36], sa[108]]], [[sa[12], sa[84]], [sa[48], sa[120]]], [[sa[24],
                                                                                           sa[96]], [sa[60], sa[132]]]],
           [[[sa[4], sa[76]], [sa[40], sa[112]]], [[sa[16],
                                                    sa[88]], [sa[52], sa[124]]],
            [[sa[28], sa[100]], [sa[64], sa[136]]]],
           [[[sa[8],
              sa[80]], [sa[44], sa[116]]], [[sa[20], sa[92]], [sa[56], sa[128]]], [[sa[32],
                                                                                    sa[104]], [sa[68], sa[140]]]]],
          [[[[sa[2], sa[74]], [sa[38], sa[110]]], [[sa[14],
                                                    sa[86]], [sa[50], sa[122]]], [[sa[26], sa[98]], [sa[62], sa[134]]]],
           [[[sa[6],
              sa[78]], [sa[42], sa[114]]], [[sa[18], sa[90]], [sa[54], sa[126]]], [[sa[30],
                                                                                    sa[102]], [sa[66], sa[138]]]],
           [[[sa[10], sa[82]], [sa[46], sa[118]]], [[sa[22],
                                                     sa[94]], [sa[58], sa[130]]],
            [[sa[34], sa[106]], [sa[70], sa[142]]]]]],
         [[[[[sa[1],
              sa[73]], [sa[37], sa[109]]], [[sa[13], sa[85]], [sa[49], sa[121]]], [[sa[25],
                                                                                    sa[97]], [sa[61], sa[133]]]],
           [[[sa[5], sa[77]], [sa[41], sa[113]]], [[sa[17],
                                                    sa[89]], [sa[53], sa[125]]],
            [[sa[29], sa[101]], [sa[65], sa[137]]]],
           [[[sa[9],
              sa[81]], [sa[45], sa[117]]], [[sa[21], sa[93]], [sa[57], sa[129]]], [[sa[33],
                                                                                    sa[105]], [sa[69], sa[141]]]]],
          [[[[sa[3], sa[75]], [sa[39], sa[111]]], [[sa[15],
                                                    sa[87]], [sa[51], sa[123]]], [[sa[27], sa[99]], [sa[63], sa[135]]]],
           [[[sa[7],
              sa[79]], [sa[43], sa[115]]], [[sa[19], sa[91]], [sa[55], sa[127]]], [[sa[31],
                                                                                    sa[103]], [sa[67], sa[139]]]],
           [[[sa[11], sa[83]], [sa[47], sa[119]]], [[sa[23],
                                                     sa[95]], [sa[59], sa[131]]],
            [[sa[35], sa[107]], [sa[71], sa[143]]]]]]])

    assert permutedims(po, (1, 0, 2, 3, 4, 5)) == Array(
        [[[[[[sa[0], sa[1]], [sa[2], sa[3]]], [[sa[4], sa[5]], [sa[6], sa[7]]], [[sa[8], sa[9]], [sa[10],
                                                                                                  sa[11]]]],
           [[[sa[12], sa[13]], [sa[14], sa[15]]], [[sa[16], sa[17]], [sa[18],
                                                                      sa[19]]], [[sa[20], sa[21]], [sa[22], sa[23]]]],
           [[[sa[24], sa[25]], [sa[26],
                                sa[27]]], [[sa[28], sa[29]], [sa[30], sa[31]]], [[sa[32], sa[33]], [sa[34],
                                                                                                    sa[35]]]]],
          [[[[sa[72], sa[73]], [sa[74], sa[75]]], [[sa[76], sa[77]], [sa[78],
                                                                      sa[79]]], [[sa[80], sa[81]], [sa[82], sa[83]]]],
           [[[sa[84], sa[85]], [sa[86],
                                sa[87]]], [[sa[88], sa[89]], [sa[90], sa[91]]], [[sa[92], sa[93]], [sa[94],
                                                                                                    sa[95]]]],
           [[[sa[96], sa[97]], [sa[98], sa[99]]], [[sa[100], sa[101]], [sa[102],
                                                                        sa[103]]],
            [[sa[104], sa[105]], [sa[106], sa[107]]]]]], [[[[[sa[36], sa[37]], [sa[38],
                                                                                sa[39]]],
                                                            [[sa[40], sa[41]], [sa[42], sa[43]]],
                                                            [[sa[44], sa[45]], [sa[46],
                                                                                sa[47]]]],
                                                           [[[sa[48], sa[49]], [sa[50], sa[51]]],
                                                            [[sa[52], sa[53]], [sa[54],
                                                                                sa[55]]],
                                                            [[sa[56], sa[57]], [sa[58], sa[59]]]],
                                                           [[[sa[60], sa[61]], [sa[62],
                                                                                sa[63]]],
                                                            [[sa[64], sa[65]], [sa[66], sa[67]]],
                                                            [[sa[68], sa[69]], [sa[70],
                                                                                sa[71]]]]], [
                                                              [[[sa[108], sa[109]], [sa[110], sa[111]]],
                                                               [[sa[112], sa[113]], [sa[114],
                                                                                     sa[115]]],
                                                               [[sa[116], sa[117]], [sa[118], sa[119]]]],
                                                              [[[sa[120], sa[121]], [sa[122],
                                                                                     sa[123]]],
                                                               [[sa[124], sa[125]], [sa[126], sa[127]]],
                                                               [[sa[128], sa[129]], [sa[130],
                                                                                     sa[131]]]],
                                                              [[[sa[132], sa[133]], [sa[134], sa[135]]],
                                                               [[sa[136], sa[137]], [sa[138],
                                                                                     sa[139]]],
                                                               [[sa[140], sa[141]], [sa[142], sa[143]]]]]]])

    assert permutedims(po, (0, 2, 1, 4, 3, 5)) == Array(
        [[[[[[sa[0], sa[1]], [sa[4], sa[5]], [sa[8], sa[9]]], [[sa[2], sa[3]], [sa[6], sa[7]], [sa[10],
                                                                                                sa[11]]]],
           [[[sa[36], sa[37]], [sa[40], sa[41]], [sa[44], sa[45]]], [[sa[38],
                                                                      sa[39]], [sa[42], sa[43]], [sa[46], sa[47]]]]],
          [[[[sa[12], sa[13]], [sa[16],
                                sa[17]], [sa[20], sa[21]]], [[sa[14], sa[15]], [sa[18], sa[19]], [sa[22],
                                                                                                  sa[23]]]],
           [[[sa[48], sa[49]], [sa[52], sa[53]], [sa[56], sa[57]]], [[sa[50],
                                                                      sa[51]], [sa[54], sa[55]], [sa[58], sa[59]]]]],
          [[[[sa[24], sa[25]], [sa[28],
                                sa[29]], [sa[32], sa[33]]], [[sa[26], sa[27]], [sa[30], sa[31]], [sa[34],
                                                                                                  sa[35]]]],
           [[[sa[60], sa[61]], [sa[64], sa[65]], [sa[68], sa[69]]], [[sa[62],
                                                                      sa[63]], [sa[66], sa[67]], [sa[70], sa[71]]]]]],
         [[[[[sa[72], sa[73]], [sa[76],
                                sa[77]], [sa[80], sa[81]]], [[sa[74], sa[75]], [sa[78], sa[79]], [sa[82],
                                                                                                  sa[83]]]],
           [[[sa[108], sa[109]], [sa[112], sa[113]], [sa[116], sa[117]]], [[sa[110],
                                                                            sa[111]], [sa[114], sa[115]],
                                                                           [sa[118], sa[119]]]]],
          [[[[sa[84], sa[85]], [sa[88],
                                sa[89]], [sa[92], sa[93]]], [[sa[86], sa[87]], [sa[90], sa[91]], [sa[94],
                                                                                                  sa[95]]]],
           [[[sa[120], sa[121]], [sa[124], sa[125]], [sa[128], sa[129]]], [[sa[122],
                                                                            sa[123]], [sa[126], sa[127]],
                                                                           [sa[130], sa[131]]]]],
          [[[[sa[96], sa[97]], [sa[100],
                                sa[101]], [sa[104], sa[105]]], [[sa[98], sa[99]], [sa[102], sa[103]], [sa[106],
                                                                                                       sa[107]]]],
           [[[sa[132], sa[133]], [sa[136], sa[137]], [sa[140], sa[141]]], [[sa[134],
                                                                            sa[135]], [sa[138], sa[139]],
                                                                           [sa[142], sa[143]]]]]]])
