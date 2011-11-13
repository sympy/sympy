""" rewrite of lambdify - This stuff is not stable at all.

It is for internal use in the new plotting module.
It may (will! see the Q'n'A in the source) be rewritten."""

import re
from sympy.printing.lambdarepr import lambdarepr
from sympy import Symbol, NumberSymbol, I, zoo, oo

#  We parse the expression string into a tree that identifies functions. Then
# we translate the names of the functions and we translate also some strings
# that are not names of functions (all this according to translation
# dictionaries).
#  If the translation goes to another module (like numpy) the
# module is imported and 'func' is translated to 'module.func'.
#  If a function can not be translated, the inner nodes of that part of the
# tree are not translated. So if we have Integral(sqrt(x)), sqrt is not
# translated to np.sqrt and the Integral does not crash.

#  Please, if there is a bug, do not try to fix it here! Rewrite this by using
# the method proposed in the last Q'n'A below.
#  If you insist on fixing it here, look at the workarounds in the function
# sympy_expression_namespace and in lambdify.

# Q: Why are you not using python abstract syntax tree?
# A: Because it is more complicated and not much more powerful in this case.

# Q: What if I have Symbol('sin') or g=Function('f')?
# A: You will break the algorithm. We should use srepr to defend against this?
#  The problem with Symbol('sin') is that it will be printed as 'sin'. The
# parser will distinguish it from the function 'sin' because functions are
# detected thanks to the opening parenthesis, but the lambda expression won't
# understand the difference if we have also the sin function.
# The solution (complicated) is to use srepr and maybe ast.
#  The problem with the g=Function('f') is that it will be printed as 'f' but in
# the global namespace we have only 'g'. But as the same printer is used in the
# constructor of the namespace there will be no problem.

# Q: What if some of the printers are not printing as expected?
# A: The algorithm wont work. You must use srepr for those cases. But even
# srepr may not print well. All problems with printers should be considered
# bugs.

# Q: What about _imp_ functions?
# A: Those are taken care for by evalf.

# Q: Will ast fix all possible problems?
# A: No. You will always have to use some printer. Even srepr may not work in
# some cases. But if the printer does not work, that should be considered a
# bug.

# Q: Is there same way to fix all possible problems?
# A: Probably by constructing our strings ourself by traversing the (func,
# args) tree and creating the namespace at the same time. That actually sounds
# good.

def experimental_lambdify(args, expr):
    # Constructing the argument string
    if not all([isinstance(a, Symbol) for a in args]):
        raise ValueError('The arguments must be Symbols.')
    else:
        argstr = ', '.join([lambdarepr(a) for a in args])


    # Constructing the translation dictionaries and making the translation
    dict_tuple_str = get_dict_tuple_str()
    dict_tuple_fun = get_dict_tuple_fun()
    exprstr = lambdarepr(expr)
    newexpr = tree2str_translate(str2tree(exprstr), dict_tuple_str, dict_tuple_fun)

    # Constructing the namespaces
    namespace = {}
    namespace.update(sympy_atoms_namespace(expr))
    namespace.update(sympy_expression_namespace(expr))
    # Ugly workaround because Pow(a,Half) prints as sqrt(a)
    # and sympy_expression_namespace can not catch it.
    from sympy import sqrt
    namespace.update({'sqrt': sqrt})
    # End workaround.
    try:
        namespace.update({'np' : __import__('numpy')})
    except ImportError:
        raise ImportError('the new plotting module needs numpy (actually the text backend does not)')


    print newexpr
    return eval('lambda ' + argstr + ' : (' + newexpr + ')', namespace)


##############################################################################
# Dicts for translating from sympy to other modules
##############################################################################

###
# numpy
###

# Functions that are the same in numpy
numpy_functions_same = [
        'sin', 'cos', 'tan',
        'sinh', 'cosh', 'tanh',
        'floor',
        'conjugate',
        ]

# Functions with different names in numpy
numpy_functions_different = {
        "acos":"arccos",
        "acosh":"arccosh",
        "arg":"angle",
        "asin":"arcsin",
        "asinh":"arcsinh",
        "atan":"arctan",
        "atan2":"arctan2",
        "atanh":"arctanh",
        "ceiling":"ceil",
        "im":"imag",
        "ln":"log",
        "Max":"amax",
        "Min":"amin",
        "re":"real",
        }

# Strings that should be translated
numpy_not_functions = {
        'I' : '1j',
        'oo':'inf',
        'E':'e',
        }

###
# Create the final ordered tuples of dictionaries
###

# For strings
def get_dict_tuple_str():
    return (numpy_not_functions,)

# For functions
def get_dict_tuple_fun():
    numpy_dict = {}
    for s in numpy_functions_same:
        numpy_dict[s] = 'np.'+s
    for k in numpy_functions_different:
        numpy_dict[k] = 'np.'+numpy_functions_different[k]
    return (numpy_dict, )

##############################################################################
# The translator functions, tree parsers, etc.
##############################################################################

def str2tree(exprstr):
    """Converts an expression string to a tree.

    Functions are represented by ('func_name(', tree_of_arguments).
    Other expressions are (head_string, mid_tree, tail_str).
    Expressions that do not contain functions are directly returned."""
    #matches the first 'function_name('
    first_par = re.match(r'[^\(]*?[\W]?(\w+\()', exprstr)
    if first_par is None:
        return exprstr
    else:
        start = len(first_par.group()) - len(first_par.groups()[0])
        head = exprstr[:start]
        func = exprstr[start:first_par.end()]
        tail = exprstr[first_par.end():]
        count = 0
        for i, c in enumerate(tail):
            if c == '(':
                count += 1
            elif c == ')':
                count -= 1
            if count == -1:
                break
        func_tail = str2tree(tail[:i])
        tail = str2tree(tail[i:])
        return (head, (func, func_tail), tail)

def tree2str(tree):
    """Converts a tree to string without translations."""
    if isinstance(tree, str):
        return tree
    else:
        return ''.join(map(tree2str, tree))

def tree2str_translate(tree, dict_tuple_str, dict_tuple_fun):
    """Converts a tree to string with translations.

    Function names are translated by translate_func.
    Other strings are translated by translate_str.
    """
    if isinstance(tree, str):
        return translate_str(tree, dict_tuple_str)
    elif isinstance(tree, tuple) and len(tree) == 2:
        return translate_func(tree[0][:-1], tree[1], dict_tuple_str, dict_tuple_fun)
    else:
        return ''.join([tree2str_translate(t, dict_tuple_str, dict_tuple_fun) for t in tree])

def translate_str(estr, dict_tuple_str):
    """Translate substrings of estr using in order the dictionaries in
    dict_tuple_str."""
    for trans_dict in dict_tuple_str:
        for k in trans_dict.keys():
            while estr.find(k) is not -1:
                i = estr.find(k)
                estr = estr[:i] + trans_dict[k] + estr[i+len(k):]
    return estr

def translate_func(func_name, argtree, dict_tuple_str, dict_tuple_fun):
    """Translate function names and the tree of arguments.

    If the function name is not in the dictionaries of dict_tuple_fun then the
    function is surrounded by an (...).evalf()."""
    for trans_dict in dict_tuple_fun:
        if func_name in trans_dict:
            new_name = trans_dict[func_name]
            argstr = tree2str_translate(argtree, dict_tuple_str, dict_tuple_fun)
            break
    else:
        new_name = '(' + func_name
        argstr = tree2str(argtree) + ').evalf()'
    return new_name + '(' + argstr

##############################################################################
# The namespace constructors
##############################################################################

def sympy_expression_namespace(expr):
    if expr is None:
        return {}
        print 'here'
    else:
        funcname = lambdarepr(expr.func)
        # Here we add an ugly workaround because lambdarepr(func(x))
        # is not always the same as lambdarepr(func). Eg
        # >>> lambdarepr(Integral(x))
        # "Integral(x)"
        # >>> lambdarepr(Integral)
        # "<class 'sympy.integrals.integrals.Integral'>"
        # >>> lambdarepr(sqrt(x))
        # "sqrt(x)"
        # >>> lambdarepr(sqrt)
        # "<function sqrt at 0x3d92de8>"
        # >>> lambdarepr(sin(x))
        # "sin(x)"
        # >>> lambdarepr(sin)
        # "sin"
        # Either one of those can be used but not all at the same time.
        # The code considers the sin example as the right one.
        regexlist = [
                r'<class \'sympy[\w.]*?.([\w]*)\'>$', # the example Integral
                r'<function ([\w]*) at 0x[\w]*>$',    # the example sqrt
                ]
        for r in regexlist:
            m = re.match(r, funcname)
            if m is not None:
                funcname = m.groups()[0]
        print funcname
        # End of the workaround
        args_dict = {}
        for a in expr.args:
            if (isinstance(a, Symbol) or
                isinstance(a, NumberSymbol) or
                a is I or a is zoo):
                continue
            else:
                args_dict.update(sympy_expression_namespace(a))
        args_dict.update({funcname : expr.func})
        return args_dict

def sympy_atoms_namespace(expr):
    atoms = expr.atoms(Symbol, NumberSymbol, I, zoo, oo)
    d = {}
    for a in atoms:
        print 'atom:' + lambdarepr(a)
        d[lambdarepr(a)] = a
    return d
