"""Printing subsystem driver"""

class Printer(object):
    """Generic printer driver

    This is a generic printer driver.
    It's job is to provide infrastructure for implementing new printers easily.

    Basically, if you want to implement a printer, all you have to do is:

    1. Subclass Printer.
    2. In your subclass, define _print_<CLASS> methods

       For each class you want to provide printing to, define an appropriate
       method how to do it. For example if you want a class FOO to be printed in
       it's own way, define _print_FOO:

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

    3. If convenient, override self.emptyPrinter

       This callable will be called to obtain printing result as a last resort,
       that is when no appropriate _print_<CLASS> was found for an expression.

    """
    def __init__(self):
        self._depth = -1
        self._str = str
        self.emptyPrinter = str

    def doprint(self, expr):
        """Returns printer's representation for expr (as a string)"""
        return self._str(self._print(expr))

    def _print(self, expr, *args):
        """internal dispatcher

           It's job is to loop through expr classes (class + it's bases), and
           try to dispatch the work to _print_<EXPR_CLASS>

           e.g., suppose we have the following class hierarchy::

                 Basic
                   |
                 Atom
                   |
                 Number
                   |
                Rational

           then, for expr=Rational(...), in order to dispatch, we will try
           calling printer methods as shown in the figure below:

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

           if nothing exists, we just return:

               p.emptyPrinter(expr)
        """
        self._depth += 1

        # See if the class of expr is known, or if one of its super
        # classes is known, and use that print function
        res = None
        for cls in type(expr).__mro__:
            if hasattr(self, '_print_'+cls.__name__):
                res = getattr(self, '_print_'+cls.__name__)(expr, *args)
                break

        # Unknown object, just use its string representation
        if res is None:
            res = self.emptyPrinter(expr)

        self._depth -= 1
        return res

