"""Printing subsystem driver

SymPy's printing system works the following way: Any expression can be
passed to a designated Printer who then is responsible to return an
adequate representation of that expression.

**The basic concept is the following:**
  1. Let the object print itself if it knows how.
  2. Take the best fitting method defined in the printer.
  3. As fall-back use the emptyPrinter method for the printer.

Which Method is Responsible for Printing?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The whole printing process is started by calling ``.doprint(expr)`` on the printer
which you want to use. This method looks for an appropriate method which can
print the given expression in the given style that the printer defines.
While looking for the method, it follows these steps:

1. **Let the object print itself if it knows how.**

    The printers looks for a specific method in every object. The name of that method
    depends on the specific printer and is defined under ``Printer.printmethod``.
    For example, StrPrinter calls ``_sympystr`` and LatexPrinter calls ``_latex``.
    Look at the documentation of the printer that you want to use.
    The name of the method is specified there.

    This was the original way of doing printing in sympy. Every class had
    its own latex, mathml, str and repr methods, but it turned out that it
    is hard to produce a high quality printer, if all the methods are spread
    out that far. Therefore all printing code was combined into the different
    printers, which works great for built-in sympy objects, but not that
    good for user defined classes where it is inconvenient to patch the
    printers.

2. **Take the best fitting method defined in the printer.**

    The printer loops through expr classes (class + its bases), and tries
    to dispatch the work to ``_print_<EXPR_CLASS>``

    e.g., suppose we have the following class hierarchy::

            Basic
            |
            Atom
            |
            Number
            |
        Rational

    then, for ``expr=Rational(...)``, the Printer will try
    to call printer methods in the order as shown in the figure below::

        p._print(expr)
        |
        |-- p._print_Rational(expr)
        |
        |-- p._print_Number(expr)
        |
        |-- p._print_Atom(expr)
        |
        `-- p._print_Basic(expr)

    if ``._print_Rational`` method exists in the printer, then it is called,
    and the result is returned back. Otherwise, the printer tries to call
    ``._print_Number`` and so on.

3. **As a fall-back use the emptyPrinter method for the printer.**

    As fall-back ``self.emptyPrinter`` will be called with the expression. If
    not defined in the Printer subclass this will be the same as ``str(expr)``.

Example of Custom Printer and Custom Printing Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _printer_example:

This example shows how you can define a custom Printer and a custom printing method
for your custom SymPy object. See the code below::

    from sympy import Basic, Function, Symbol
    from sympy.printing.printer import Printer
    from sympy.printing.latex import print_latex
    from sympy.core.basic import Basic

    class MyBasic(Basic):
        \"\"\" Our custom SymPy object. You can also subclass sympy.core.Expr
        or sympy.core.Function. Both of them are subclasses of Basic.
        \"\"\"

        def __init__(self, text):
            self.text = text

        def __str__(self):
            return self.text

        def _latex(self, printer=None):
            # This is the printmethod of LatexPrinter.
            return \"\\\\textbf{{{}}}\".format(self.text)

    class MyPrinter(Printer):
        \"\"\" Our custom Printer. \"\"\"

        # This method is called for every SymPy object when we call doprint()
        # on the printer.
        printmethod = '_myprinter'

        def _print_Derivative(self, expr):
            # expr.args == (x(t), t)
            # expr.args[0] == x(t)
            # expr.args[0].func == x
            return str(expr.args[0].func) + \"'\" * len(expr.args[1:])

        def _print_MyBasic(self, expr):
            # Since MyBasic does not implement _myprinter method, this method
            # will be called instead.
            return str(expr).upper()

    def my_printer(expr):
        \"\"\" Most of the printers define their own wrappers for print().
        These wrappers usually take printer settings. Our printer does not have
        any settings.
        \"\"\"
        print(MyPrinter().doprint(expr))

    t = Symbol('t')
    x = Function('x')(t)
    dxdt = x.diff(t)
    d2xdt2 = dxdt.diff(t)

    ex = MyBasic('My expression on steroids.')

    print(d2xdt2)
    my_printer(d2xdt2)

    print(dxdt)
    my_printer(dxdt)

    print(ex)
    my_printer(ex)

    # This is the wrapper of print() for LatexPrinter.
    print_latex(ex)

The output of the code above is::

    Derivative(x(t), t, t)
    x''
    Derivative(x(t), t)
    x'
    My expression on steroids.
    MY EXPRESSION ON STEROIDS.
    \\textbf{My expression on steroids.}

"""

from __future__ import print_function, division

from sympy import Basic, Add

from sympy.core.core import BasicMeta

from functools import cmp_to_key


class Printer(object):
    """ Generic printer

    Its job is to provide infrastructure for implementing new printers easily.

    If you want to define your custom Printer or your custom printing method
    for your custom class then see the example above: printer_example_ .
    """

    _global_settings = {}

    _default_settings = {}

    emptyPrinter = str
    printmethod = None

    def __init__(self, settings=None):
        self._str = str

        self._settings = self._default_settings.copy()

        for key, val in self._global_settings.items():
            if key in self._default_settings:
                self._settings[key] = val

        if settings is not None:
            self._settings.update(settings)

            if len(self._settings) > len(self._default_settings):
                for key in self._settings:
                    if key not in self._default_settings:
                        raise TypeError("Unknown setting '%s'." % key)

        # _print_level is the number of times self._print() was recursively
        # called. See StrPrinter._print_Float() for an example of usage
        self._print_level = 0

    @classmethod
    def set_global_settings(cls, **settings):
        """Set system-wide printing settings. """
        for key, val in settings.items():
            if val is not None:
                cls._global_settings[key] = val

    @property
    def order(self):
        if 'order' in self._settings:
            return self._settings['order']
        else:
            raise AttributeError("No order defined.")

    def doprint(self, expr):
        """Returns printer's representation for expr (as a string)"""
        return self._str(self._print(expr))

    def _print(self, expr, *args, **kwargs):
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
                return getattr(expr, self.printmethod)(self, *args, **kwargs)

            # See if the class of expr is known, or if one of its super
            # classes is known, and use that print function
            for cls in type(expr).__mro__:
                printmethod = '_print_' + cls.__name__
                if hasattr(self, printmethod):
                    return getattr(self, printmethod)(expr, *args, **kwargs)
            # Unknown object, fall back to the emptyPrinter.
            return self.emptyPrinter(expr)
        finally:
            self._print_level -= 1

    def _as_ordered_terms(self, expr, order=None):
        """A compatibility function for ordering terms in Add. """
        order = order or self.order

        if order == 'old':
            return sorted(Add.make_args(expr), key=cmp_to_key(Basic._compare_pretty))
        else:
            return expr.as_ordered_terms(order=order)
