"""Printing subsystem driver

SymPy's printing system works the following way: Any expression can be passed to
a designated Printer who then is responsible to return a adequate representation
of that expression.

The basic concept is the following:
  1. Let the object print itself if it knows how.
  2. Take the best fitting method defined in the printer.
  3. As fall-back use the emptyPrinter method for the printer.

Some more information how the single concepts work and who should use which:

1. The object prints itself

    This was the original way of doing printing in sympy. Every class had its
    own latex, mathml, str and repr methods, but it turned out that it is hard
    to produce a high quality printer, if all the methods are spread out that
    far. Therefor all printing code was combined into the different printers,
    which works great for built-in sympy objects, but not that good for user
    defined classes where it is inconvenient to patch the printers.
    To get nevertheless a fitting representation, the printers look for a
    specific method in every object, that will be called if it's available and
    is then responsible for the representation. The name of that method depends
    on the specific printer and is defined under Printer.printmethodname.


2. Take the best fitting method defined in the printer.

    The printer loops through expr classes (class + its bases), and tries to dispatch the
    work to _print_<EXPR_CLASS>

    e.g., suppose we have the following class hierarchy::

            Basic
            |
            Atom
            |
            Number
            |
        Rational

    then, for expr=Rational(...), in order to dispatch, we will try calling printer methods
    as shown in the figure below::

        p._print(expr)
        |
        |-- p._print_Rational(expr)
        |
        |-- p._print_Number(expr)
        |
        |-- p._print_Atom(expr)
        |
        `-- p._print_Basic(expr)

    if ._print_Rational method exists in the printer, then it is called,
    and the result is returned back.

    otherwise, we proceed with trying Rational bases in the inheritance
    order.

3. As fall-back use the emptyPrinter method for the printer.

    As fall-back self.emptyPrinter will be called with the expression. If
    not defined in the Printer subclass this will be the same as str(expr)
"""

from sympy import Basic, Mul, Add

from sympy.core.exprtools import decompose_power
from sympy.polys.monomialtools import monomial_cmp
from sympy.core.basic import BasicMeta

class Printer(object):
    """Generic printer

    Its job is to provide infrastructure for implementing new printers easily.

    Basically, if you want to implement a printer, all you have to do is:

    1. Subclass Printer.

    2. Define Printer.printmethod in your subclass.
       If a object has a method with that name, this method will be used
       for printing.

    3. In your subclass, define _print_<CLASS> methods

       For each class you want to provide printing to, define an appropriate
       method how to do it. For example if you want a class FOO to be printed in
       its own way, define _print_FOO:

       def _print_FOO(self, e):
           ...

       this should return how FOO instance e is printed

       Also, if BAR is a subclass of FOO, _print_FOO(bar) will be called for
       instance of BAR, if no _print_BAR is provided.  Thus, usually, we don't
       need to provide printing routines for every class we want to support --
       only generic routine has to be provided for a set of classes.

       A good example for this are functions - for example PrettyPrinter only
       defines _print_Function, and there is no _print_sin, _print_tan, etc...

       On the other hand, a good printer will probably have to define separate
       routines for Symbol, Atom, Number, Integral, Limit, etc...

    4. If convenient, override self.emptyPrinter

       This callable will be called to obtain printing result as a last resort,
       that is when no appropriate print method was found for an expression.

    Example of overloading StrPrinter::

        from sympy import Basic, Function, Symbol
        from sympy.printing.str import StrPrinter

        class CustomStrPrinter(StrPrinter):
            \"\"\"
            Example of how to customize the StrPrinter for both a Sympy class and a
            user defined class subclassed from the Sympy Basic class.
            \"\"\"

            def _print_Derivative(self, expr):
                \"\"\"
                Custom printing of the Sympy Derivative class.

                Instead of:

                D(x(t), t) or D(x(t), t, t)

                We will print:

                x'     or     x''

                In this example, expr.args == (x(t), t), and expr.args[0] == x(t), and
                expr.args[0].func == x
                \"\"\"
                return str(expr.args[0].func) + "'"*len(expr.args[1:])

            def _print_MyClass(self, expr):
                \"\"\"
                Print the characters of MyClass.s alternatively lower case and upper
                case
                \"\"\"
                s = ""
                i = 0
                for char in expr.s:
                    if i % 2 == 0:
                        s += char.lower()
                    else:
                        s += char.upper()
                    i += 1
                return s

        # Override the __str__ method of to use CustromStrPrinter
        Basic.__str__ = lambda self: CustomStrPrinter().doprint(self)
        # Demonstration of CustomStrPrinter:
        t = Symbol('t')
        x = Function('x')(t)
        dxdt = x.diff(t)            # dxdt is a Derivative instance
        d2xdt2 = dxdt.diff(t)       # dxdt2 is a Derivative instance
        ex = MyClass('I like both lowercase and upper case')

        print dxdt
        print d2xdt2
        print ex

    The output of the above code is::

        x'
        x''
        i lIkE BoTh lOwErCaSe aNd uPpEr cAsE

    By overriding Basic.__str__, we can customize the printing of anything that
    is subclassed from Basic.

    """

    _default_settings = {
        "order": None,
    }

    emptyPrinter = str
    printmethod = None

    def __init__(self, settings=None):
        self._str = str
        self._settings = self._default_settings.copy()
        if settings is not None:
            self._settings.update(settings)
            if len(self._settings) > len(self._default_settings):
                for key in self._settings:
                    if key not in self._default_settings:
                        raise TypeError("Unknown setting '%s'." % key)

        # _print_level is the number of times self._print() was recursively
        # called. See StrPrinter._print_Real() for an example of usage
        self._print_level = 0

        # configure ordering of terms in Add (e.g. lexicographic ordering)
        self.order, self.reverse = self.parse_order(self._settings["order"])

    def doprint(self, expr):
        """Returns printer's representation for expr (as a string)"""
        return self._str(self._print(expr))

    def _print(self, expr, *args):
        """Internal dispatcher

        Tries the following concepts to print an expression:
            1. Let the object print itself if it knows how.
            2. Take the best fitting method defined in the printer.
            3. As fall-back use the emptyPrinter method for the printer.
        """
        self._print_level += 1
        try:
            # If the printer defines a name for a printing method
            # (Printer.printmethod) and the object knows for itself how it
            # should be printed, use that method.
            if (self.printmethod and hasattr(expr, self.printmethod)
                    and not isinstance(expr, BasicMeta)):
                return getattr(expr, self.printmethod)(self, *args)

            # See if the class of expr is known, or if one of its super
            # classes is known, and use that print function
            for cls in type(expr).__mro__:
                printmethod = '_print_' + cls.__name__
                if hasattr(self, printmethod):
                    return getattr(self, printmethod)(expr, *args)

            # Unknown object, fall back to the emptyPrinter.
            return self.emptyPrinter(expr)
        finally:
            self._print_level -= 1

    def analyze(self, expr, order=None):
        """Rewrite an expression as sorted list of terms. """
        if order is not None:
            order, reverse = self.parse_order(order)
        else:
            order, reverse = None, None

        gens, terms = set([]), []

        for term in Add.make_args(expr):
            coeff, cpart, ncpart = [], {}, []

            for factor in Mul.make_args(term):
                if not factor.is_commutative:
                    ncpart.append(factor)
                else:
                    if factor.is_Number:
                        coeff.append(factor)
                    else:
                        base, exp = decompose_power(factor)

                        if exp < 0:
                            exp, base = -exp, 1/base

                        cpart[base] = exp
                        gens.add(base)

            terms.append((coeff, cpart, ncpart, term))

        gens = sorted(gens, Basic._compare_pretty)

        k, indices = len(gens), {}

        for i, g in enumerate(gens):
            indices[g] = i

        result = []

        for coeff, cpart, ncpart, term in terms:
            monom = [0]*k

            for base, exp in cpart.iteritems():
                monom[indices[base]] = exp

            result.append((coeff, monom, ncpart, term))

        if order is None and self.order is None:
            return sorted(result, Basic._compare_pretty)
        else:
            return sorted(result, lambda a, b: \
                self._compare_terms(a, b, order=order, reverse=reverse))

    def _compare_terms(self, a, b, order=None, reverse=None):
        """Compare two terms using data from Printer.analyze(). """
        a_coeff, a_monom, a_ncpart, _ = a
        b_coeff, b_monom, b_ncpart, _ = b

        if order is None:
            order = self.order
        if reverse is None:
            reverse = self.reverse

        result = order(a_monom, b_monom)

        if not result:
            if not (a_ncpart or b_ncpart):
                result = cmp(a_coeff, b_coeff)
            else:
                result = Basic._compare_pretty(Mul(*a_ncpart), Mul(*b_ncpart))

                if not result:
                    result = cmp(a_coeff, b_coeff)

        if not reverse:
            return -result
        else:
            return result

    @classmethod
    def parse_order(cls, order):
        """Parse and configure ordering of terms in Add. """
        if order is not None:
            if order.startswith('rev-'):
                order = monomial_cmp(order[4:])
                reverse = True
            else:
                order = monomial_cmp(order)
                reverse = False
        else:
            order = None
            reverse = False

        return order, reverse
