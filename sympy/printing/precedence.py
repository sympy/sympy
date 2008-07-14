"""A module providing information about the necessity of brackets"""

# Default precedence values for some basic types
PRECEDENCE = {
    "Lambda":1,
    "Relational":20,
    "Add":40,
    "Mul":50,
    "Pow":60,
    "Apply":70,
    "Item":75,
    "Atom":1000
}

# A dictionary assigning precedence values to certain classes. These values are
# treated like they were inherited, so not every single class has to be named
# here.
PRECEDENCE_VALUES = {}

# Sometimes it's not enough to assign a fixed precedence value to a
# class. Then a function can be inserted in this dictionary that takes
# an instance of this class as argument and returns the appropriate
# precedence value.
PRECEDENCE_FUNCTIONS = {}

# Clearly arranged list of classes with their precedence level. These classes
# are imported and inserted into PRECEDENCE_VALUES.
_PREC_VALS = (
    ("sympy.core.add.Add", "Add"),
    ("sympy.core.basic.Atom", "Atom"),
    ("sympy.core.function.Function", "Apply"),
    ("sympy.core.function.Derivative", "Apply"),
    ("sympy.core.numbers.Integer", "Atom"),
    ("sympy.core.power.Pow", "Pow"),
    ("sympy.core.relational.Relational", "Relational"),
    ("sympy.simplify.cse_opts.Sub", "Add"),
    ("sympy.series.order.Order", "Apply"),
    ("sympy.integrals.integrals.Integral", "Apply"),
    ("sympy.concrete.products.Product", "Apply"),
    ("sympy.concrete.summations.Sum", "Apply"),
)

# Insert information of _PREC_VALS propperly into PRECEDENCE_VALUES
for imp, level in _PREC_VALS:
    mod, cls = imp.rsplit(".", 1)
    mod = __import__(mod, fromlist=[cls])
    PRECEDENCE_VALUES[getattr(mod, cls)] = PRECEDENCE[level]


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

# Clearly arranged list of classes with their precedence functions. These
#  functions are inserted into PRECEDENCE_FUNCTIONS.
_PREC_FUNCS =(
    ("sympy.core.mul.Mul", precedence_Mul),
    ("sympy.core.numbers.Rational", precedence_Rational),
)

# Insert information of _PREC_FUNCS propperly into PRECEDENCE_FUNCS
for imp, func in _PREC_FUNCS:
    mod, cls = imp.rsplit(".", 1)
    mod = __import__(mod, fromlist=[cls])
    PRECEDENCE_FUNCTIONS[getattr(mod, cls)] = func


def precedence(item):
    """
    Returns the precedence of a given object.
    """
    if hasattr(item, "precedence"):
            return item.precedence
    for i in item.__class__.__mro__:
        if i in PRECEDENCE_FUNCTIONS:
            return PRECEDENCE_FUNCTIONS[i](item)
        elif i in PRECEDENCE_VALUES:
            return PRECEDENCE_VALUES[i]
    return 1000
