import sympy.physics.mechanics.eombase as eombase

# Determine what inputs should be used
EOM = eombase.EOMBase(frame=None)

# Should contain access to all the property attributes contained in existing EOM
# generators

Thing = "Place holder until a correct output is chosen"

# assert EOM.auxiliary_eqs == Thing
assert EOM.bodies == Thing
assert EOM.mass_matrix == Thing
assert EOM.mass_matrix_full == Thing
assert EOM.loads == Thing
assert EOM.forcing == Thing
assert EOM.forcing_full == Thing
assert EOM.coordinates == Thing
assert EOM.speeds == Thing

# The properties should not be able to be altered
# EOM.auxiliary_eqs = Thing  # Raises error
EOM.bodies = Thing  # Raises error
EOM.mass_matrix = Thing  # Raises error
EOM.mass_matrix_full = Thing  # Raises error
EOM.loads = Thing  # Raises error
EOM.forcing = Thing  # Raises error
EOM.forcing_full = Thing  # Raises error
EOM.coordinates = Thing  # Raises error
EOM.speeds = Thing  # Raises error
