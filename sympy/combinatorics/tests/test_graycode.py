from sympy.combinatorics.graycode import GrayCode

def test_graycode():
    a = GrayCode(6)
    assert len(list(a.doit())) == 64
    assert list(a.doit('011001')) == ['011001', '011011', '011010',
    '011110', '011111', '011101', '011100', '010100', '010101', '010111',
    '010110', '010010', '010011', '010001', '010000', '110000', '110001',
    '110011', '110010', '110110', '110111', '110101', '110100', '111100',
    '111101', '111111', '111110', '111010', '111011', '111001', '111000',
    '101000', '101001', '101011', '101010', '101110', '101111', '101101',
    '101100', '100100', '100101', '100111', '100110', '100010', '100011',
    '100001', '100000']

    a = GrayCode(5, start = '10010')
    assert a.rank_current == 28
    a = GrayCode(6, start = '101000')
    assert a.rank_current == 48

    assert GrayCode.unrank_gray(4, 6).current == '000110'
    assert GrayCode.unrank_gray(4, 6).rank_current == 4
    assert [GrayCode(4, start=s).rank_current for s in \
            GrayCode(4).doit()] == [0, 1, 2, 3, 4, 5, 6, 7, 8, \
                                    9, 10, 11, 12, 13, 14, 15]
    a = GrayCode(15, rank=15)
    assert a.current == '000000000001000'

