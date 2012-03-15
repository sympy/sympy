from sympy import inversej,sni,I


def test_inversej():
    assert inversej('sn',1.5,0.3) == 1.71388944817879 + 1.27343861737965*I
    assert inversej('sn',0.5,0.7) == 0.540872222535946
    assert inversej('ns',0.5,0.7) == 0.672614860709727 + 1.71388944817879*I
    assert inversej('dn',0.5,0.7) == 2.07536313529247 + 0.496697108138456*I
    assert inversej('nd',0.5,0.7) == 1.18325254863892*I
    assert inversej('cn',0.5,0.7) == 1.19663065156446
    assert inversej('nc',0.5,0.7) == 1.09913522309204*I
    assert inversej('sc',0.5,0.7) == 0.475562686531585
    assert inversej('cs',0.5,0.7) == 1.28530363049307
    assert inversej('cd',0.5,0.7) == -1.53449091275652
    assert inversej('dc',0.5,0.7) == -1.40274827458274 + 1.71388944817879*I
    assert inversej('sd',0.5,0.7) == 3.65811036106799
    assert inversej('ds',0.5,0.7) == 2.07536313529247 + 0.424207693683713*I
def test_sni():
    assert sni(4+3*I,0.6)   == 0.204791309225623 + 1.61832120477254*I
    assert sni(8+.3*I,0.6)  == 0.162273290847638 + 1.77134810630362*I
    assert sni(4.4*I,0.6)   == 1.49057438078435*I
    assert sni(5+4.4*I,0.6) == 0.144394206350047 + 1.64787545453375*I
