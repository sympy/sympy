from sympy.liealgebras.cartan_type import CartanType

def test_type_E():
    c = CartanType("E6")
    assert c.dimension() == 6
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "2\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for i in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(3, 7))
    assert c.dynkin_diagram() == diag
    posroots = c.positive_roots()
    assert posroots[8] == [1, 0, 0, 0, 1, 0, 0, 0]
