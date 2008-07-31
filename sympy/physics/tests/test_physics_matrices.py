from sympy.physics.matrices import msigma, mgamma, minkowski_tensor
from sympy import zeros, eye, I



def test_Pauli():
    #this and the following test are testing both Pauli and Dirac matrices
    #and also that the general Matrix class works correctly in a real world
    #situation
    sigma1=msigma(1)
    sigma2=msigma(2)
    sigma3=msigma(3)

    assert sigma1 == sigma1
    assert sigma1 != sigma2

    # sigma*I -> I*sigma    (see #354)
    assert sigma1*sigma2 == sigma3*I
    assert sigma3*sigma1 == sigma2*I
    assert sigma2*sigma3 == sigma1*I

    assert sigma1*sigma1 == eye(2)
    assert sigma2*sigma2 == eye(2)
    assert sigma3*sigma3 == eye(2)

    assert sigma1*2*sigma1 == 2*eye(2)
    assert sigma1*sigma3*sigma1 == -sigma3

def test_Dirac():
    gamma0=mgamma(0)
    gamma1=mgamma(1)
    gamma2=mgamma(2)
    gamma3=mgamma(3)
    gamma5=mgamma(5)

    # gamma*I -> I*gamma    (see #354)
    assert gamma5 == gamma0 * gamma1 * gamma2 * gamma3 * I
    assert gamma1 * gamma2 + gamma2 * gamma1 == zeros(4)
    assert gamma0 * gamma0 == eye(4) * minkowski_tensor[0,0]
    assert gamma2 * gamma2 != eye(4) * minkowski_tensor[0,0]
    assert gamma2 * gamma2 == eye(4) * minkowski_tensor[2,2]

    assert mgamma(5,True) == \
        mgamma(0,True)*mgamma(1,True)*mgamma(2,True)*mgamma(3,True)*I


