from sympy.physics.mechanics import MovingRefFrame
####################################################################################

#The base frame for all frame trees in the tests
O = MovingRefFrame('O')

def test_pos_vector():
    """
    Tests for relative position vectors of reference frames
    """
    
