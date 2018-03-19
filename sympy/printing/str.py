"""
A Printer for generating readable representation of most sympy classes.
"""

from __future__ import print_function, division

from sympy.core import S, Rational, Pow, Basic, Mul
from sympy.core.mul import _keep_coeff
from .printer import Printer
from sympy.printing.precedence import precedence, PRECEDENCE

import mpmath.libmp as mlib
from mpmath.libmp import prec_to_dps

from sympy.utilities import default_sort_key


class StrPrinter(Printer):
    printmethod = "_sympystr"
    _default_settings = {
        "order": None,
        "full_prec": "auto",
        "sympy_integers": False,
        "abbrev": False,
    }

    _relationals = dict()

    def parenthesize(self, item, level, strict=False, **kwargs):
        if (precedence(item) < level) or ((not strict) and precedence(item) <= level):
            return "(%s)" % self._print(item, **kwargs)
        else:
            return self._print(item, **kwargs)

    def stringify(self, args, sep, level=0):
        return sep.join([self.parenthesize(item, level) for item in args])

    def emptyPrinter(self, expr):
        if isinstance(expr, str):
            return expr
        elif isinstance(expr, Basic):
            if hasattr(expr, "args"):
                return repr(expr)
            else:
                raise
        else:
            return str(expr)

    def _print_Add(self, expr, **kwargs):
        order = kwargs.get('order', None)
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = self._as_ordered_terms(expr, order=order)

        PREC = precedence(expr)
        l = []
        for term in terms:
            t = self._print(term, **kwargs)
            if t.startswith('-'):
                sign = "-"
                t = t[1:]
            else:
                sign = "+"
            if precedence(term) < PREC:
                l.extend([sign, "(%s)" % t])
            else:
                l.extend([sign, t])
        sign = l.pop(0)
        if sign == '+':
            sign = ""
        return sign + ' '.join(l)

    def _print_BooleanTrue(self, expr, **kwargs):
        return "True"

    def _print_BooleanFalse(self, expr, **kwargs):
        return "False"

    def _print_Not(self, expr, **kwargs):
        return '~%s' %(self.parenthesize(expr.args[0],PRECEDENCE["Not"]))

    def _print_And(self, expr, **kwargs):
        return self.stringify(expr.args, " & ", PRECEDENCE["BitwiseAnd"])

    def _print_Or(self, expr, **kwargs):
        return self.stringify(expr.args, " | ", PRECEDENCE["BitwiseOr"])

    def _print_AppliedPredicate(self, expr, **kwargs):
        return '%s(%s)' % (self._print(expr.func, **kwargs), self._print(expr.arg, **kwargs))

    def _print_Basic(self, expr, **kwargs):
        l = [self._print(o, **kwargs) for o in expr.args]
        return expr.__class__.__name__ + "(%s)" % ", ".join(l)

    def _print_BlockMatrix(self, B, **kwargs):
        if B.blocks.shape == (1, 1):
            self._print(B.blocks[0, 0], **kwargs)
        return self._print(B.blocks, **kwargs)

    def _print_Catalan(self, expr, **kwargs):
        return 'Catalan'

    def _print_ComplexInfinity(self, expr, **kwargs):
        return 'zoo'

    def _print_Derivative(self, expr, **kwargs):
        dexpr = expr.expr
        dvars = [i[0] if i[1] == 1 else i for i in expr.variable_count]
        return 'Derivative(%s)' % ", ".join(map(lambda arg: self._print(arg, **kwargs), [dexpr] + dvars))

    def _print_dict(self, d, **kwargs):
        keys = sorted(d.keys(), key=default_sort_key)
        items = []

        for key in keys:
            item = "%s: %s" % (self._print(key, **kwargs), self._print(d[key], **kwargs))
            items.append(item)

        return "{%s}" % ", ".join(items)

    def _print_Dict(self, expr, **kwargs):
        return self._print_dict(expr, **kwargs)

    def _print_RandomDomain(self, d, **kwargs):
        if hasattr(d, 'as_boolean'):
            return 'Domain: ' + self._print(d.as_boolean(), **kwargs)
        elif hasattr(d, 'set'):
            return ('Domain: ' + self._print(d.symbols, **kwargs) + ' in ' +
                    self._print(d.set, **kwargs))
        else:
            return 'Domain on ' + self._print(d.symbols)

    def _print_Dummy(self, expr, **kwargs):
        return '_' + expr.name

    def _print_EulerGamma(self, expr, **kwargs):
        return 'EulerGamma'

    def _print_Exp1(self, expr, **kwargs):
        return 'E'

    def _print_ExprCondPair(self, expr, **kwargs):
        return '(%s, %s)' % (self._print(expr.expr, **kwargs), self._print(expr.cond, **kwargs))

    def _print_FiniteSet(self, s, **kwargs):
        s = sorted(s, key=default_sort_key)
        if len(s) > 10:
            printset = s[:3] + ['...'] + s[-3:]
        else:
            printset = s
        return '{' + ', '.join(self._print(el, **kwargs) for el in printset) + '}'

    def _print_Function(self, expr, **kwargs):
        return expr.func.__name__ + "(%s)" % self.stringify(expr.args, ", ")

    def _print_GeometryEntity(self, expr, **kwargs):
        # GeometryEntity is special -- it's base is tuple
        return str(expr)

    def _print_GoldenRatio(self, expr, **kwargs):
        return 'GoldenRatio'

    def _print_ImaginaryUnit(self, expr, **kwargs):
        return 'I'

    def _print_Infinity(self, expr, **kwargs):
        return 'oo'

    def _print_Integral(self, expr, **kwargs):
        def _xab_tostr(xab):
            if len(xab) == 1:
                return self._print(xab[0], **kwargs)
            else:
                return self._print((xab[0],) + tuple(xab[1:]), **kwargs)
        L = ', '.join([_xab_tostr(l) for l in expr.limits])
        return 'Integral(%s, %s)' % (self._print(expr.function, **kwargs), L)

    def _print_Interval(self, i, **kwargs):
        fin =  'Interval{m}({a}, {b})'
        a, b, l, r = i.args
        if a.is_infinite and b.is_infinite:
            m = ''
        elif a.is_infinite and not r:
            m = ''
        elif b.is_infinite and not l:
            m = ''
        elif not l and not r:
            m = ''
        elif l and r:
            m = '.open'
        elif l:
            m = '.Lopen'
        else:
            m = '.Ropen'
        return fin.format(**{'a': a, 'b': b, 'm': m})

    def _print_AccumulationBounds(self, i, **kwargs):
        return "AccumBounds(%s, %s)" % (self._print(i.min, **kwargs),
                                        self._print(i.max, **kwargs))

    def _print_Inverse(self, I, **kwargs):
        return "%s^-1" % self.parenthesize(I.arg, PRECEDENCE["Pow"])

    def _print_Lambda(self, obj, **kwargs):
        args, expr = obj.args
        if len(args) == 1:
            return "Lambda(%s, %s)" % (self._print(args.args[0]), self._print(expr))
        else:
            arg_string = ", ".join(self._print(arg, **kwargs) for arg in args)
            return "Lambda((%s), %s)" % (arg_string, self._print(expr, **kwargs))

    def _print_LatticeOp(self, expr, **kwargs):
        args = sorted(expr.args, key=default_sort_key)
        return expr.func.__name__ + "(%s)" % ", ".join(self._print(arg, **kwargs) for arg in args)

    def _print_Limit(self, expr, **kwargs):
        e, z, z0, dir = expr.args
        if str(dir) == "+":
            return "Limit(%s, %s, %s)" % tuple(map(self._print, (e, z, z0)))
        else:
            return "Limit(%s, %s, %s, dir='%s')" % tuple(map(self._print,
                                                            (e, z, z0, dir)))

    def _print_list(self, expr, **kwargs):
        return "[%s]" % self.stringify(expr, ", ")

    def _print_MatrixBase(self, expr, **kwargs):
        return expr._format_str(self)
    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_MatrixElement(self, expr, **kwargs):
        return self.parenthesize(expr.parent, PRECEDENCE["Atom"], strict=True) \
            + '[%s, %s]' % (self._print(expr.i), self._print(expr.j))

    def _print_MatrixSlice(self, expr, **kwargs):
        def strslice(x):
            x = list(x)
            if x[2] == 1:
                del x[2]
            if x[1] == x[0] + 1:
                del x[1]
            if x[0] == 0:
                x[0] = ''
            return ':'.join(map(lambda arg: self._print(arg, **kwargs), x))
        return (self._print(expr.parent, **kwargs) + '[' +
                strslice(expr.rowslice) + ', ' +
                strslice(expr.colslice) + ']')

    def _print_DeferredVector(self, expr, **kwargs):
        return expr.name

    def _print_Mul(self, expr, **kwargs):

        prec = precedence(expr)

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append(Rational(item.p))
                if item.q != 1:
                    b.append(Rational(item.q))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec, strict=False) for x in a]
        b_str = [self.parenthesize(x, prec, strict=False) for x in b]

        if len(b) == 0:
            return sign + '*'.join(a_str)
        elif len(b) == 1:
            return sign + '*'.join(a_str) + "/" + b_str[0]
        else:
            return sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str)

    def _print_MatMul(self, expr, **kwargs):
        c, m = expr.as_coeff_mmul()
        if c.is_number and c < 0:
            expr = _keep_coeff(-c, m)
            sign = "-"
        else:
            sign = ""

        return sign + '*'.join(
            [self.parenthesize(arg, precedence(expr)) for arg in expr.args]
        )

    def _print_HadamardProduct(self, expr, **kwargs):
        return '.*'.join([self.parenthesize(arg, precedence(expr))
            for arg in expr.args])

    def _print_MatAdd(self, expr, **kwargs):
        terms = [self.parenthesize(arg, precedence(expr))
             for arg in expr.args]
        l = []
        for t in terms:
            if t.startswith('-'):
                sign = "-"
                t = t[1:]
            else:
                sign = "+"
            l.extend([sign, t])
        sign = l.pop(0)
        if sign == '+':
            sign = ""
        return sign + ' '.join(l)

    def _print_NaN(self, expr, **kwargs):
        return 'nan'

    def _print_NegativeInfinity(self, expr, **kwargs):
        return '-oo'

    def _print_Normal(self, expr, **kwargs):
        return "Normal(%s, %s)" % (self._print(expr.mu, **kwargs), self._print(expr.sigma, **kwargs))

    def _print_Order(self, expr, **kwargs):
        if all(p is S.Zero for p in expr.point) or not len(expr.variables):
            if len(expr.variables) <= 1:
                return 'O(%s)' % self._print(expr.expr, **kwargs)
            else:
                return 'O(%s)' % self.stringify((expr.expr,) + expr.variables, ', ', 0)
        else:
            return 'O(%s)' % self.stringify(expr.args, ', ', 0)

    def _print_Ordinal(self, expr, **kwargs):
        return expr.__str__()

    def _print_Cycle(self, expr, **kwargs):
        return expr.__str__()

    def _print_Permutation(self, expr, **kwargs):
        from sympy.combinatorics.permutations import Permutation, Cycle
        if Permutation.print_cyclic:
            if not expr.size:
                return '()'
            # before taking Cycle notation, see if the last element is
            # a singleton and move it to the head of the string
            s = Cycle(expr)(expr.size - 1).__repr__()[len('Cycle'):]
            last = s.rfind('(')
            if not last == 0 and ',' not in s[last:]:
                s = s[last:] + s[:last]
            s = s.replace(',', '')
            return s
        else:
            s = expr.support()
            if not s:
                if expr.size < 5:
                    return 'Permutation(%s)' % self._print(expr.array_form)
                return 'Permutation([], size=%s)' % self._print(expr.size)
            trim = self._print(expr.array_form[:s[-1] + 1]) + ', size=%s' % self._print(expr.size)
            use = full = self._print(expr.array_form)
            if len(trim) < len(full):
                use = trim
            return 'Permutation(%s)' % use

    def _print_TensorIndex(self, expr, **kwargs):
        return expr._print(**kwargs)

    def _print_TensorHead(self, expr, **kwargs):
        return expr._print(**kwargs)

    def _print_Tensor(self, expr, **kwargs):
        return expr._print(**kwargs)

    def _print_TensMul(self, expr, **kwargs):
        return expr._print(**kwargs)

    def _print_TensAdd(self, expr, **kwargs):
        return expr._print(**kwargs)

    def _print_PermutationGroup(self, expr, **kwargs):
        p = ['    %s' % self._print(a, **kwargs) for a in expr.args]
        return 'PermutationGroup([\n%s])' % ',\n'.join(p)

    def _print_PDF(self, expr, **kwargs):
        return 'PDF(%s, (%s, %s, %s))' % \
            (self._print(expr.pdf.args[1], **kwargs), self._print(expr.pdf.args[0], **kwargs),
            self._print(expr.domain[0], **kwargs), self._print(expr.domain[1], **kwargs))

    def _print_Pi(self, expr, **kwargs):
        return 'pi'

    def _print_PolyRing(self, ring, **kwargs):
        return "Polynomial ring in %s over %s with %s order" % \
            (", ".join(map(lambda rs: self._print(rs, **kwargs), ring.symbols)),
            self._print(ring.domain, **kwargs), self._print(ring.order, **kwargs))

    def _print_FracField(self, field, **kwargs):
        return "Rational function field in %s over %s with %s order" % \
            (", ".join(map(lambda fs: self._print(fs, **kwargs), field.symbols)),
            self._print(field.domain, **kwargs), self._print(field.order, **kwargs))

    def _print_FreeGroupElement(self, elm, **kwargs):
        return elm.__str__()

    def _print_PolyElement(self, poly, **kwargs):
        return poly.str(self, PRECEDENCE, "%s**%s", "*")

    def _print_FracElement(self, frac, **kwargs):
        if frac.denom == 1:
            return self._print(frac.numer, **kwargs)
        else:
            numer = self.parenthesize(frac.numer, PRECEDENCE["Mul"], strict=True)
            denom = self.parenthesize(frac.denom, PRECEDENCE["Atom"], strict=True)
            return numer + "/" + denom

    def _print_Poly(self, expr, **kwargs):
        ATOM_PREC = PRECEDENCE["Atom"] - 1
        terms, gens = [], [ self.parenthesize(s, ATOM_PREC) for s in expr.gens ]

        for monom, coeff in expr.terms():
            s_monom = []

            for i, exp in enumerate(monom):
                if exp > 0:
                    if exp == 1:
                        s_monom.append(gens[i])
                    else:
                        s_monom.append(gens[i] + "**%d" % exp)

            s_monom = "*".join(s_monom)

            if coeff.is_Add:
                if s_monom:
                    s_coeff = "(" + self._print(coeff, **kwargs) + ")"
                else:
                    s_coeff = self._print(coeff, **kwargs)
            else:
                if s_monom:
                    if coeff is S.One:
                        terms.extend(['+', s_monom])
                        continue

                    if coeff is S.NegativeOne:
                        terms.extend(['-', s_monom])
                        continue

                s_coeff = self._print(coeff, **kwargs)

            if not s_monom:
                s_term = s_coeff
            else:
                s_term = s_coeff + "*" + s_monom

            if s_term.startswith('-'):
                terms.extend(['-', s_term[1:]])
            else:
                terms.extend(['+', s_term])

        if terms[0] in ['-', '+']:
            modifier = terms.pop(0)

            if modifier == '-':
                terms[0] = '-' + terms[0]

        format = expr.__class__.__name__ + "(%s, %s"

        from sympy.polys.polyerrors import PolynomialError

        try:
            format += ", modulus=%s" % expr.get_modulus()
        except PolynomialError:
            format += ", domain='%s'" % expr.get_domain()

        format += ")"

        for index, item in enumerate(gens):
            if len(item) > 2 and (item[:1] == "(" and item[len(item) - 1:] == ")"):
                gens[index] = item[1:len(item) - 1]

        return format % (' '.join(terms), ', '.join(gens))

    def _print_ProductSet(self, p, **kwargs):
        return ' x '.join(self._print(set, **kwargs) for set in p.sets)

    def _print_AlgebraicNumber(self, expr, **kwargs):
        if expr.is_aliased:
            return self._print(expr.as_poly().as_expr(), **kwargs)
        else:
            return self._print(expr.as_expr(), **kwargs)

    def _print_Pow(self, expr, **kwargs):
        rational = kwargs.get('rational', False)
        PREC = precedence(expr)

        if expr.exp is S.Half and not rational:
            return "sqrt(%s)" % self._print(expr.base, **kwargs)

        if expr.is_commutative:
            if -expr.exp is S.Half and not rational:
                # Note: Don't test "expr.exp == -S.Half" here, because that will
                # match -0.5, which we don't want.
                return "%s/sqrt(%s)" % tuple(map(lambda arg: self._print(arg, **kwargs), (S.One, expr.base)))
            if expr.exp is -S.One:
                # Similarly to the S.Half case, don't test with "==" here.
                return '%s/%s' % (self._print(S.One, **kwargs),
                                  self.parenthesize(expr.base, PREC, strict=False))

        e = self.parenthesize(expr.exp, PREC, strict=False)
        if self.printmethod == '_sympyrepr' and expr.exp.is_Rational and expr.exp.q != 1:
            # the parenthesized exp should be '(Rational(a, b))' so strip parens,
            # but just check to be sure.
            if e.startswith('(Rational'):
                return '%s**%s' % (self.parenthesize(expr.base, PREC, strict=False), e[1:-1])
        return '%s**%s' % (self.parenthesize(expr.base, PREC, strict=False), e)

    def _print_UnevaluatedExpr(self, expr, **kwargs):
        return self._print(expr.args[0], **kwargs)

    def _print_MatPow(self, expr, **kwargs):
        PREC = precedence(expr)
        return '%s**%s' % (self.parenthesize(expr.base, PREC, strict=False),
                         self.parenthesize(expr.exp, PREC, strict=False))

    def _print_ImmutableDenseNDimArray(self, expr, **kwargs):
        return str(expr)

    def _print_ImmutableSparseNDimArray(self, expr, **kwargs):
        return str(expr)

    def _print_Integer(self, expr, **kwargs):
        if self._settings.get("sympy_integers", False):
            return "S(%s)" % (expr)
        return str(expr.p)

    def _print_Integers(self, expr, **kwargs):
        return 'S.Integers'

    def _print_Naturals(self, expr, **kwargs):
        return 'S.Naturals'

    def _print_Naturals0(self, expr, **kwargs):
        return 'S.Naturals0'

    def _print_Reals(self, expr, **kwargs):
        return 'S.Reals'

    def _print_int(self, expr, **kwargs):
        return str(expr)

    def _print_mpz(self, expr, **kwargs):
        return str(expr)

    def _print_Rational(self, expr, **kwargs):
        if expr.q == 1:
            return str(expr.p)
        else:
            if self._settings.get("sympy_integers", False):
                return "S(%s)/%s" % (expr.p, expr.q)
            return "%s/%s" % (expr.p, expr.q)

    def _print_PythonRational(self, expr, **kwargs):
        if expr.q == 1:
            return str(expr.p)
        else:
            return "%d/%d" % (expr.p, expr.q)

    def _print_Fraction(self, expr, **kwargs):
        if expr.denominator == 1:
            return str(expr.numerator)
        else:
            return "%s/%s" % (expr.numerator, expr.denominator)

    def _print_mpq(self, expr, **kwargs):
        if expr.denominator == 1:
            return str(expr.numerator)
        else:
            return "%s/%s" % (expr.numerator, expr.denominator)

    def _print_Float(self, expr, **kwargs):
        prec = expr._prec
        if prec < 5:
            dps = 0
        else:
            dps = prec_to_dps(expr._prec)
        if self._settings["full_prec"] is True:
            strip = False
        elif self._settings["full_prec"] is False:
            strip = True
        elif self._settings["full_prec"] == "auto":
            strip = self._print_level > 1
        rv = mlib.to_str(expr._mpf_, dps, strip_zeros=strip)
        if rv.startswith('-.0'):
            rv = '-0.' + rv[3:]
        elif rv.startswith('.0'):
            rv = '0.' + rv[2:]
        if rv.startswith('+'):
            # e.g., +inf -> inf
            rv = rv[1:]
        return rv

    def _print_Relational(self, expr, **kwargs):

        charmap = {
            "==": "Eq",
            "!=": "Ne",
            ":=": "Assignment",
            '+=': "AddAugmentedAssignment",
            "-=": "SubAugmentedAssignment",
            "*=": "MulAugmentedAssignment",
            "/=": "DivAugmentedAssignment",
            "%=": "ModAugmentedAssignment",
        }

        if expr.rel_op in charmap:
            return '%s(%s, %s)' % (charmap[expr.rel_op], self._print(expr.lhs),
                                   self._print(expr.rhs))

        return '%s %s %s' % (self.parenthesize(expr.lhs, precedence(expr)),
                           self._relationals.get(expr.rel_op) or expr.rel_op,
                           self.parenthesize(expr.rhs, precedence(expr)))

    def _print_ComplexRootOf(self, expr, **kwargs):
        return "CRootOf(%s, %d)" % (self._print_Add(expr.expr,  order='lex', **kwargs),
                                    expr.index)

    def _print_RootSum(self, expr, **kwargs):
        args = [self._print_Add(expr.expr, order='lex', **kwargs)]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun, **kwargs))

        return "RootSum(%s)" % ", ".join(args)

    def _print_GroebnerBasis(self, basis, **kwargs):
        cls = basis.__class__.__name__

        exprs = [self._print_Add(arg, order=basis.order, **kwargs) for arg in basis.exprs]
        exprs = "[%s]" % ", ".join(exprs)

        gens = [ self._print(gen, **kwargs) for gen in basis.gens ]
        domain = "domain='%s'" % self._print(basis.domain, **kwargs)
        order = "order='%s'" % self._print(basis.order, **kwargs)

        args = [exprs] + gens + [domain, order]

        return "%s(%s)" % (cls, ", ".join(args))

    def _print_Sample(self, expr, **kwargs):
        return "Sample([%s])" % self.stringify(expr, ", ", 0)

    def _print_set(self, s, **kwargs):
        items = sorted(s, key=default_sort_key)

        args = ', '.join(self._print(item, **kwargs) for item in items)
        if not args:
            return "set()"
        return '{%s}' % args

    def _print_frozenset(self, s, **kwargs):
        if not s:
            return "frozenset()"
        return "frozenset(%s)" % self._print_set(s, **kwargs)

    def _print_SparseMatrix(self, expr, **kwargs):
        from sympy.matrices import Matrix
        return self._print(Matrix(expr), **kwargs)

    def _print_Sum(self, expr, **kwargs):
        def _xab_tostr(xab):
            if len(xab) == 1:
                return self._print(xab[0], **kwargs)
            else:
                return self._print((xab[0],) + tuple(xab[1:]), **kwargs)
        L = ', '.join([_xab_tostr(l) for l in expr.limits])
        return 'Sum(%s, %s)' % (self._print(expr.function, **kwargs), L)

    def _print_Symbol(self, expr, **kwargs):
        return expr.name
    _print_MatrixSymbol = _print_Symbol
    _print_RandomSymbol = _print_Symbol

    def _print_Identity(self, expr, **kwargs):
        return "I"

    def _print_ZeroMatrix(self, expr, **kwargs):
        return "0"

    def _print_Predicate(self, expr, **kwargs):
        return "Q.%s" % expr.name

    def _print_str(self, expr, **kwargs):
        return str(expr)

    def _print_tuple(self, expr, **kwargs):
        if len(expr) == 1:
            return "(%s,)" % self._print(expr[0], **kwargs)
        else:
            return "(%s)" % self.stringify(expr, ", ")

    def _print_Tuple(self, expr, **kwargs):
        return self._print_tuple(expr, **kwargs)

    def _print_Transpose(self, T, **kwargs):
        return "%s.T" % self.parenthesize(T.arg, PRECEDENCE["Pow"])

    def _print_Uniform(self, expr, **kwargs):
        return "Uniform(%s, %s)" % (self._print(expr.a, **kwargs), self._print(expr.b, **kwargs))

    def _print_Union(self, expr, **kwargs):
        return 'Union(%s)' %(', '.join([self._print(a, **kwargs) for a in expr.args]))

    def _print_Complement(self, expr, **kwargs):
        return r' \ '.join(self._print(set_, **kwargs) for set_ in expr.args)

    def _print_Quantity(self, expr, **kwargs):
        if self._settings.get("abbrev", False):
            return "%s" % expr.abbrev
        return "%s" % expr.name

    def _print_Quaternion(self, expr, **kwargs):
        s = [self.parenthesize(i, PRECEDENCE["Mul"], strict=True) for i in expr.args]
        a = [s[0]] + [i+"*"+j for i, j in zip(s[1:], "ijk")]
        return " + ".join(a)

    def _print_Dimension(self, expr, **kwargs):
        return str(expr)

    def _print_Wild(self, expr, **kwargs):
        return expr.name + '_'

    def _print_WildFunction(self, expr, **kwargs):
        return expr.name + '_'

    def _print_Zero(self, expr, **kwargs):
        if self._settings.get("sympy_integers", False):
            return "S(0)"
        return "0"

    def _print_DMP(self, p, **kwargs):
        from sympy.core.sympify import SympifyError
        try:
            if p.ring is not None:
                # TODO incorporate order
                return self._print(p.ring.to_sympy(p), **kwargs)
        except SympifyError:
            pass

        cls = p.__class__.__name__
        rep = self._print(p.rep, **kwargs)
        dom = self._print(p.dom, **kwargs)
        ring = self._print(p.ring, **kwargs)

        return "%s(%s, %s, %s)" % (cls, rep, dom, ring)

    def _print_DMF(self, expr, **kwargs):
        return self._print_DMP(expr, **kwargs)

    def _print_Object(self, obj, **kwargs):
        return 'Object("%s")' % obj.name

    def _print_IdentityMorphism(self, morphism, **kwargs):
        return 'IdentityMorphism(%s)' % morphism.domain

    def _print_NamedMorphism(self, morphism, **kwargs):
        return 'NamedMorphism(%s, %s, "%s")' % \
               (morphism.domain, morphism.codomain, morphism.name)

    def _print_Category(self, category, **kwargs):
        return 'Category("%s")' % category.name

    def _print_BaseScalarField(self, field, **kwargs):
        return field._coord_sys._names[field._index]

    def _print_BaseVectorField(self, field, **kwargs):
        return 'e_%s' % field._coord_sys._names[field._index]

    def _print_Differential(self, diff, **kwargs):
        field = diff._form_field
        if hasattr(field, '_coord_sys'):
            return 'd%s' % field._coord_sys._names[field._index]
        else:
            return 'd(%s)' % self._print(field, **kwargs)

    def _print_Tr(self, expr, **kwargs):
        #TODO : Handle indices
        return "%s(%s)" % ("Tr", self._print(expr.args[0], **kwargs))


def sstr(expr, **settings):
    """Returns the expression as a string.

    For large expressions where speed is a concern, use the setting
    order='none'. If abbrev=True setting is used then units are printed in
    abbreviated form.

    Examples
    ========

    >>> from sympy import symbols, Eq, sstr
    >>> a, b = symbols('a b')
    >>> sstr(Eq(a + b, 0))
    'Eq(a + b, 0)'
    """

    p = StrPrinter(settings)
    s = p.doprint(expr)

    return s


class StrReprPrinter(StrPrinter):
    """(internal) -- see sstrrepr"""

    def _print_str(self, s):
        return repr(s)


def sstrrepr(expr, **settings):
    """return expr in mixed str/repr form

       i.e. strings are returned in repr form with quotes, and everything else
       is returned in str form.

       This function could be useful for hooking into sys.displayhook
    """

    p = StrReprPrinter(settings)
    s = p.doprint(expr)

    return s
