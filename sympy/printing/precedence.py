"""A module providing information about the necessity of brackets"""

# Default precedence values for some basic types
PRECEDENCE = {
    "Lambda":1,
    "Relational":20,
    "Or":20,
    "And":30,
    "Add":40,
    "Mul":50,
    "Pow":60,
    "Not":100,
    "Atom":1000
}

# A dictionary assigning precedence values to certain classes. These values are
# treated like they were inherited, so not every single class has to be named
# here.
PRECEDENCE_VALUES = {
    "Or" : PRECEDENCE["Or"],
    "And" : PRECEDENCE["And"],
    "Add" : PRECEDENCE["Add"],
    "Pow" : PRECEDENCE["Pow"],
    "Relational" : PRECEDENCE["Relational"],
    "Sub" : PRECEDENCE["Add"],
    "Not": PRECEDENCE["Not"],
}

# Sometimes it's not enough to assign a fixed precedence value to a
# class. Then a function can be inserted in this dictionary that takes
# an instance of this class as argument and returns the appropriate
# precedence value.

# Precedence functions
def precedence_Mul(item):
    coeff, rest = item.as_coeff_terms()
    if coeff.is_negative:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Mul"]

def precedence_Rational(item):
    if item.p < 0:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Mul"]

def precedence_Integer(item):
    if item.p < 0:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Atom"]

PRECEDENCE_FUNCTIONS = {
    "Integer" : precedence_Integer,
    "Mul" : precedence_Mul,
    "Rational" : precedence_Rational,
}


def precedence(item):
    """
    Returns the precedence of a given object.
    """
    if hasattr(item, "precedence"):
        return item.precedence
    for i in item.__class__.__mro__:
        n = i.__name__
        if n in PRECEDENCE_FUNCTIONS:
            return PRECEDENCE_FUNCTIONS[n](item)
        elif n in PRECEDENCE_VALUES:
            return PRECEDENCE_VALUES[n]
    return PRECEDENCE["Atom"]
