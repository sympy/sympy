# Fix imports to work in both direct run and python -m
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from sympy.physics.units.aviation_units import parse_flight_level, convert_quantity_to_numeric
from sympy.physics.units import meter, foot

def test_flight_level():
    # Test valid FL
    fl = parse_flight_level("FL290")
    print("FL290 in meters:", convert_quantity_to_numeric(fl, meter))  # ~8839.2
    print("FL290 in feet:", convert_quantity_to_numeric(fl, foot))     # 29000

    # Test out-of-range
    try:
        parse_flight_level("FL2000")
    except ValueError as e:
        print("Error:", e)

if __name__ == "__main__":
    test_flight_level()
