from sympy.combinatorics.graycode import GrayCode

def test_graycode():
    a = GrayCode(4)
    assert len(list(a.generate_bitlist())) == 16
    assert list(a.generate_bitlist(['0', '0', '1', '0'])) == \
    [['0', '0', '1', '0'], ['0', '1', '1', '0'], ['0', '1', '1', '1'], \
    ['0', '1', '0', '1'], ['0', '1', '0', '0'], ['1', '1', '0', '0'], \
    ['1', '1', '0', '1'], ['1', '1', '1', '1'], ['1', '1', '1', '0'], \
    ['1', '0', '1', '0'], ['1', '0', '1', '1'], ['1', '0', '0', '1'], \
    ['1', '0', '0', '0']]

    a = GrayCode(5, start = ['1','0','0','1','0'])
    assert a.rank == 14
    a = GrayCode(6, start = ['1','0','1','0','0','0'])
    assert a.rank == 6

    assert GrayCode.unrank_gray_code(4, 6)._current == \
          ['0', '1', '1', '0', '0', '0']
    assert GrayCode.unrank_gray_code(4, 6).rank == 4
    assert [GrayCode(4, start=s).rank for s in \
            GrayCode(4).generate_bitlist()] == [0, 15, 8, 7, 4, 11,
                                                12, 3, 2, 13, 10, 5,
                                                6, 9, 14, 1]

