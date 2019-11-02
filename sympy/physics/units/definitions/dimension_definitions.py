from sympy.physics.units import Dimension


angle = Dimension(name="angle")

# base dimensions (MKS)
length = Dimension(name="length", symbol="L")
mass = Dimension(name="mass", symbol="M")
time = Dimension(name="time", symbol="T")

# base dimensions (MKSA not in MKS)
current = Dimension(name='current', symbol='I')

# other base dimensions:
temperature = Dimension("temperature", "T")
amount_of_substance = Dimension("amount_of_substance")
luminous_intensity = Dimension("luminous_intensity")

# derived dimensions (MKS)
velocity = Dimension(name="velocity")
acceleration = Dimension(name="acceleration")
momentum = Dimension(name="momentum")
force = Dimension(name="force", symbol="F")
energy = Dimension(name="energy", symbol="E")
power = Dimension(name="power")
pressure = Dimension(name="pressure")
frequency = Dimension(name="frequency", symbol="f")
action = Dimension(name="action", symbol="A")
volume = Dimension("volume")

# derived dimensions (MKSA not in MKS)
voltage = Dimension(name='voltage', symbol='U')
impedance = Dimension(name='impedance', symbol='Z')
conductance = Dimension(name='conductance', symbol='G')
capacitance = Dimension(name='capacitance')
inductance = Dimension(name='inductance')
charge = Dimension(name='charge', symbol='Q')
magnetic_density = Dimension(name='magnetic_density', symbol='B')
magnetic_flux = Dimension(name='magnetic_flux')

# Dimensions in information theory:
information = Dimension(name='information')
