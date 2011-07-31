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
