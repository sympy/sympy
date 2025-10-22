import re
from sympy.physics.units import Quantity, foot, meter

# Flight Level unit: 1 FL = 100 ft
FL = Quantity('flight_level', factor=100 * foot)

MIN_FL = -10
MAX_FL = 1267

def parse_flight_level(fl_str):
    """
    Parse 'FLxxx' string to SymPy Quantity (flight_level)
    """
    m = re.match(r"FL(\d{1,4})", fl_str.upper())
    if not m:
        raise ValueError(f"Invalid Flight Level string: {fl_str}")

    fl_number = int(m.group(1))
    if not (MIN_FL <= fl_number <= MAX_FL):
        raise ValueError(f"Flight Level {fl_number} out of range ({MIN_FL}..{MAX_FL})")

    return FL * fl_number  # Quantity object

def convert_quantity_to_numeric(qty, unit):
    """
    Convert Quantity to numeric value in given unit (meter or foot)
    """
    return (qty / unit).evalf()
