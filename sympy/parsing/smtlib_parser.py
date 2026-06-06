"""
This module contains the classes and functions used to parse SMT-LIB 2.7
code into native SymPy expressions.
"""

from sympy.core.symbol import Symbol
from sympy.core.function import Function
from sympy.core.numbers import Integer, Float
from sympy.logic.boolalg import (
    And, Or, Not, Implies, Xor, Equivalent
)
from sympy.core.relational import (
    Eq, Ne, StrictLessThan, StrictGreaterThan, LessThan, GreaterThan
)
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.complexes import Abs


class SMTLibSyntaxError(Exception):
    """
    Exception raised for syntax errors in SMT-LIB scripts.
    """

    def __init__(self, message, line, col):
        """
        Initialize the syntax error with line and column information.

        Parameters
        ==========

        message : str
            A string containing the error message
        line : int
            An integer denoting the line number of the error
        col : int
            An integer denoting the column number of the error
        """
        super().__init__(f"{message} at line {line}, col {col}")


class UnknownSMTLibCommandError(Exception):
    """
    Exception raised for unsupported SMT-LIB commands.
    """
    pass


class SMTTokenizer:
    """
    A lexical analyzer for SMT-LIB 2.7 source code.
    """

    def __init__(self, source):
        """
        Initialize the tokenizer with the SMT-LIB source code string.

        Parameters
        ==========

        source : str
            A string containing the SMT-LIB source code
        """
        self.source = source

    def tokenize(self):
        """
        Generate lexical tokens from the source string.

        Returns
        =======

        result : iterator
            A generator yielding tuples of (token, line_number, column)
        """
        line_num = 1
        col_num = 0
        token = []
        token_start_col = 0
        # Tracks if lexer is inside a ';' line comment
        in_comment = False
        # Tracks if lexer is inside a "string literal"
        in_string = False
        in_quoted_symbol = False  # Tracks if lexer is inside a |quoted symbol|

        for char in self.source:
            if char == '\n':
                if not in_string and not in_quoted_symbol:
                    if token:
                        yield "".join(token), line_num, token_start_col
                        token = []
                    in_comment = False
                else:
                    token.append(char)
                line_num += 1
                col_num = 0
                continue

            col_num += 1
            if in_comment:
                continue

            if char == '"' and not in_quoted_symbol:
                # SMT-LIB string literals preserve spaces and semicolons.
                # We toggle state here to prevent premature token splitting.
                if not token and not in_string:
                    token_start_col = col_num
                token.append(char)
                in_string = not in_string
                if not in_string:
                    yield "".join(token), line_num, token_start_col
                    token = []
                continue

            if char == '|' and not in_string:
                if not token and not in_quoted_symbol:
                    token_start_col = col_num
                token.append(char)
                in_quoted_symbol = not in_quoted_symbol
                if not in_quoted_symbol:
                    yield "".join(token), line_num, token_start_col
                    token = []
                continue

            if in_string or in_quoted_symbol:
                token.append(char)
                continue

            if char == ';':
                if token:
                    yield "".join(token), line_num, token_start_col
                    token = []
                in_comment = True
                continue

            if char in ' \t\r':
                if token:
                    yield "".join(token), line_num, token_start_col
                    token = []
                continue

            if char in '()':
                if token:
                    yield "".join(token), line_num, token_start_col
                    token = []
                yield char, line_num, col_num
                continue

            if not token:
                token_start_col = col_num
            token.append(char)

        if token:
            yield "".join(token), line_num, token_start_col


class SMTLibParser:
    """
    Parser that converts SMT-LIB 2.7 expressions to SymPy objects.
    """

    def __init__(self):
        """
        Initialize parser state, symbol tables, and stack contexts.
        """
        # Global registry mapping SMT-LIB variable names to SymPy Symbols
        self._symbols = {}
        # Registry mapping uninterpreted function names to SymPy Function
        # objects
        self._functions = {}
        # Stores macro definitions mapping name -> (arg_names, temp_syms,
        # parsed_body)
        self._macros = {}
        # Accumulates the list of mathematical assertions generated by the
        # script
        self._assertions = []
        # Tracks the declared underlying logic structure of the SMT script
        # (e.g. QF_LRA)
        self._logic = None

        # Tracks declared custom sorts (types) mapping name -> arity
        self._sorts = {}
        # Stores generic global parser/solver configurations (e.g. :print-
        # success)
        self._options = {}
        # Stores metadata about the SMT-LIB benchmark (e.g. :author, :status)
        self._info = {}

        # Pre-registered set of operators corresponding to extended theories
        # (e.g., Strings, BitVectors). Operations hitting this set are
        # preserved
        # as generic uninterpreted functions rather than failing.
        self._extended_ops = {
            'to_real', 'to_int', 'is_int',
            'concat', 'bvnot', 'bvand', 'bvor', 'bvneg', 'bvadd', 'bvmul',
            'bvudiv', 'bvurem', 'bvshl', 'bvlshr', 'bvsub', 'bvult', 'bvxor',
            'bvnand', 'bvnor', 'bvxnor', 'bvcomp', 'bvsdiv', 'bvsrem',
            'bvsmod', 'bvashr', 'bvule', 'bvugt', 'bvuge', 'bvslt', 'bvsle',
            'bvsgt', 'bvsge', 'bv2nat', 'bv2int', 'ubv_to_int',
            'ext_rotate_left', 'ext_rotate_right', 'str.len', 'str.++',
            'str.at', 'str.contains', 'str.indexof', 'str.replace',
            'str.substr', 'str.prefixof', 'str.suffixof', 'str.to.int',
            'int.to.str', 'select', 'store', 'div', 'mod'
        }

        # Each stack maintains a historical snapshot to restore previous
        # scoping environments
        self._assertion_stack = []
        self._symbol_stack = []
        self._macro_stack = []
        self._function_stack = []
        self._sort_stack = []

    def parse_str(self, source):
        """
        Parse an SMT-LIB string and return symbols and assertions.

        Parameters
        ==========

        source : str
            A string containing the SMT-LIB source code to be parsed

        Returns
        =======

        (symbols, assertions)
            A tuple containing the symbols and assertions
        """
        sexps = self._get_s_expressions(source)
        for expr in sexps:
            self.transform(expr)
        return self._symbols, self._assertions

    def _get_s_expressions(self, source):
        """
        Lex a source string and build nested list representations.

        Parameters
        ==========

        source : str
            A string containing the SMT-LIB source code

        Returns
        =======

        result : current_list
            A nested list of parsed S-expressions
        """
        tokenizer = SMTTokenizer(source)
        stack = []
        current_list = []

        for token, line, col in tokenizer.tokenize():
            if token == '(':
                stack.append(current_list)
                current_list = []
            elif token == ')':
                if not stack:
                    raise SMTLibSyntaxError("Unexpected ')'", line, col)
                parent_list = stack.pop()
                parent_list.append(current_list)
                current_list = parent_list
            else:
                current_list.append(token)

        if stack:
            raise SMTLibSyntaxError("Unexpected EOF, missing ')'", line, col)

        return current_list

    def transform(self, expr):
        """
        Recursively transform an S-expression into a SymPy AST.

        Parameters
        ==========

        expr : list
            A parsed S-expression list or atomic token

        Returns
        =======

        result : Any
            A SymPy AST object representing the expression
        """
        if not isinstance(expr, list) or len(expr) == 0:
            return self._transform_atom(expr)

        # Normalize the command name, replacing '-' with '_' to match Python
        # methods
        command = expr[0]
        method_name = str(command).replace('-', '_').lower()
        handler = getattr(self, f'transform_{method_name}', None)

        # If no specific structural handler exists, assume it's a logical
        # operator
        if handler is None:
            return self.transform_operator(expr)

        return handler(expr)

    def _transform_atom(self, atom):
        """
        Convert a scalar atom into a SymPy type or native boolean.

        Parameters
        ==========

        atom : str
            A string representing an atomic value

        Returns
        =======

        result : Any
            A native boolean or a SymPy Integer, Float, or Symbol
        """
        if atom in self._symbols:
            return self._symbols[atom]

        atom_str = str(atom).lower()
        if atom_str == 'true':
            return True
        if atom_str == 'false':
            return False

        # Check standard SMT-LIB hex and binary numeric literals
        if atom_str.startswith('#x'):
            return Integer(int(atom[2:], 16))
        if atom_str.startswith('#b'):
            return Integer(int(atom[2:], 2))

        try:
            # Attempt to parse as a real/float or an integer
            if '.' in atom:
                return Float(atom)
            return Integer(atom)
        except ValueError:
            # If not numeric or boolean, it represents a symbolic variable
            return Symbol(str(atom))

    def transform_set_logic(self, expr):
        """
        Define the underlying logic of the SMT script.

        Parameters
        ==========

        expr : list
            A list containing the command and the logic name
        """
        self._logic = expr[1]

    def transform_declare_sort(self, expr):
        """
        Register a custom sort declaration with its arity.

        Parameters
        ==========

        expr : list
            A list containing the command, sort name, and optional arity
        """
        sort_name = expr[1]
        arity = int(expr[2]) if len(expr) > 2 else 0
        self._sorts[sort_name] = arity

    def transform_declare_sort_parameter(self, expr):
        """
        Register a sort parameter declaration with zero arity.

        Parameters
        ==========

        expr : list
            A list containing the command and the sort name
        """
        sort_name = expr[1]
        self._sorts[sort_name] = 0

    def transform_define_sort(self, expr):
        """
        Define a custom sort combining existing parameters.

        Parameters
        ==========

        expr : list
            A list containing the command, name, parameters, and definition
        """
        sort_name = expr[1]
        sort_params = expr[2]
        sort_def = expr[3]
        self._sorts[sort_name] = (sort_params, sort_def)

    def transform_check_sat(self, expr):
        """
        Process check-sat natively without invoking a solver here.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        pass

    def transform_check_sat_assuming(self, expr):
        """
        Process check-sat-assuming without invoking a solver here.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        pass

    def transform_exit(self, expr):
        """
        Gracefully process the exit command.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        pass

    def transform_set_option(self, expr):
        """
        Set a global parser or solver option.

        Parameters
        ==========

        expr : list
            A list containing the command, option name, and value
        """
        self._options[expr[1]] = expr[2] if len(expr) > 2 else True

    def transform_set_info(self, expr):
        """
        Set a global parser or solver info metadata.

        Parameters
        ==========

        expr : list
            A list containing the command, info key, and value
        """
        self._info[expr[1]] = expr[2] if len(expr) > 2 else True

    def transform_echo(self, expr):
        """
        Print the given string to standard output.

        Parameters
        ==========

        expr : str
            A list containing the command and the string to print
        """
        if len(expr) > 1 and isinstance(expr[1], str):
            print(expr[1].strip('"'))

    def transform_reset(self, expr):
        """
        Completely reset the internal parser state.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        self._symbols.clear()
        self._functions.clear()
        self._macros.clear()
        self._assertions.clear()
        self._sorts.clear()
        self._options.clear()
        self._info.clear()
        self._assertion_stack.clear()
        self._symbol_stack.clear()
        self._macro_stack.clear()
        self._function_stack.clear()
        self._sort_stack.clear()
        self._logic = None

    def transform_reset_assertions(self, expr):
        """
        Clear all assertions without affecting global declarations.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        self._assertions.clear()
        self._assertion_stack.clear()
        self._symbol_stack.clear()
        self._macro_stack.clear()
        self._function_stack.clear()
        self._sort_stack.clear()

    def transform_get_assertions(self, expr):
        """
        Print current set of generated SymPy assertions.

        Parameters
        ==========

        expr : list
            A list containing the command
        """
        print(self._assertions)

    def transform_get_option(self, expr):
        """
        Print the value of a defined option.

        Parameters
        ==========

        expr : list
            A list containing the command and the option name
        """
        opt = expr[1]
        print(self._options.get(opt, "unsupported"))

    def transform_get_info(self, expr):
        """
        Print the value of a defined info tag.

        Parameters
        ==========

        expr : list
            A list containing the command and the info tag
        """
        info_req = expr[1]
        print(self._info.get(info_req, "unsupported"))

    def transform_push(self, expr):
        """
        Push the current state bounds onto the increment stack.

        Parameters
        ==========

        expr : list
            A list containing the command and an optional depth
        """
        levels = int(expr[1]) if len(expr) > 1 else 1
        # Save a snapshot of the current environment state to the stack
        for _ in range(levels):
            self._assertion_stack.append(list(self._assertions))
            self._symbol_stack.append(dict(self._symbols))
            self._macro_stack.append(dict(self._macros))
            self._function_stack.append(dict(self._functions))
            self._sort_stack.append(dict(self._sorts))

    def transform_pop(self, expr):
        """
        Pop the current state bounds from the increment stack.

        Parameters
        ==========

        expr : list
            A list containing the command and an optional depth
        """
        levels = int(expr[1]) if len(expr) > 1 else 1
        # Restore the environment state from the stack snapshot
        for _ in range(levels):
            self._assertions = self._assertion_stack.pop()
            self._symbols = self._symbol_stack.pop()
            self._macros = self._macro_stack.pop()
            self._functions = self._function_stack.pop()
            self._sorts = self._sort_stack.pop()

    def _create_symbol(self, name, sort):
        """
        Generate a SymPy Symbol carrying strict SMT-LIB type traits.

        Parameters
        ==========

        name : str
            A string indicating the variable name
        sort : str
            A string indicating the variable SMT-LIB sort

        Returns
        =======

        result : Symbol
            A SymPy symbol customized with necessary assumptions
        """
        if sort == 'Real':
            return Symbol(name, real=True)
        if sort == 'Int':
            return Symbol(name, integer=True)
        if sort == 'Bool':
            return Symbol(name, boolean=True)
        return Symbol(name)

    def transform_declare_const(self, expr):
        """
        Register a standalone constant variable into the symbol table.

        Parameters
        ==========

        expr : list
            A list containing the command, variable name, and sort
        """
        name = expr[1]
        sort = expr[2]
        self._symbols[name] = self._create_symbol(name, sort)

    def transform_declare_fun(self, expr):
        """
        Register a function signature or zero-arity variable.

        Parameters
        ==========

        expr : list
            A list containing the command, name, arguments, and sort
        """
        name = expr[1]
        args = expr[2]
        sort = expr[3] if len(expr) > 3 else None
        if not args:
            self._symbols[name] = self._create_symbol(name, sort)
        else:
            self._functions[name] = Function(name)
            self._symbols[name] = self._functions[name]

    def transform_define_fun(self, expr):
        """
        Evaluate a parameterized macro and store it globally.

        Parameters
        ==========

        expr : list
            A list containing the command, name, arguments, sort, and body
        """
        name = expr[1]
        arg_decls = expr[2]
        body = expr[4]

        arg_names = [arg[0] for arg in arg_decls]
        temp_syms = {arg: Symbol(arg) for arg in arg_names}

        # Macro parameters must temporarily shadow global symbols of the same
        # name to ensure the macro body evaluates with correct local scoping.
        shadowed = {}
        for arg in arg_names:
            if arg in self._symbols:
                shadowed[arg] = self._symbols[arg]
            self._symbols[arg] = temp_syms[arg]

        parsed_body = self.transform(body)

        for arg in arg_names:
            if arg in shadowed:
                self._symbols[arg] = shadowed[arg]
            else:
                del self._symbols[arg]

        if not arg_names:
            self._symbols[name] = parsed_body
        else:
            self._macros[name] = (arg_names, temp_syms, parsed_body)

    def transform_define_const(self, expr):
        """
        Evaluate an expression block and store it as a global constant.

        Parameters
        ==========

        expr : list
            A list containing the command, constant name, sort, and value
        """
        name = expr[1]
        body = expr[3]
        self._symbols[name] = self.transform(body)

    def transform_define_fun_rec(self, expr):
        """
        Register a recursively defined function signature and body.

        Parameters
        ==========

        expr : list
            A list containing the command, name, arguments, sort, and body
        """
        self.transform_define_fun(expr)

    def transform_define_funs_rec(self, expr):
        """
        Iteratively process multiple interdependent recursive functions.

        Parameters
        ==========

        expr : list
            A list containing the command, declarations, and bodies
        """
        decls = expr[1]
        bodies = expr[2]
        for decl, body in zip(decls, bodies):
            mock_expr = ['define-fun', decl[0], decl[1], decl[2], body]
            self.transform_define_fun(mock_expr)

    def transform_declare_datatype(self, expr):
        """
        Register datatype constructor and selector signatures.

        Parameters
        ==========

        expr : list
            A list containing the command, name, and datatype declaration
        """
        dt_decl = expr[2]
        if isinstance(dt_decl, list) and dt_decl[0] == 'par':
            dt_decl = dt_decl[2]

        for constructor_decl in dt_decl:
            if isinstance(constructor_decl, list):
                c_name = constructor_decl[0]
                self._functions[c_name] = Function(c_name)
                self._symbols[c_name] = self._functions[c_name]
                for selector in constructor_decl[1:]:
                    s_name = selector[0]
                    self._functions[s_name] = Function(s_name)
                    self._symbols[s_name] = self._functions[s_name]
            else:
                c_name = constructor_decl
                self._symbols[c_name] = Symbol(c_name)

    def transform_declare_datatypes(self, expr):
        """
        Process batch parameterised datatypes simultaneously.

        Parameters
        ==========

        expr : list
            A list containing the command, sort declarations, and types
        """
        sort_decls = expr[1]
        datatype_decls = expr[2]
        for sort_decl, dt_decl in zip(sort_decls, datatype_decls):
            mock_expr = ['declare-datatype', sort_decl[0], dt_decl]
            self.transform_declare_datatype(mock_expr)

    def transform_assert(self, expr):
        """
        Evaluate an expression and append it to the assertions list.

        Parameters
        ==========

        expr : list
            A list containing the command and the expression logic
        """
        self._assertions.append(self.transform(expr[1]))

    def transform_let(self, expr):
        """
        Evaluate local variable bindings surrounding an expression block.

        Parameters
        ==========

        expr : list
            A list containing the command, local bindings, and body

        Returns
        =======

        result : Any
            The AST of the body with bounds substituted
        """
        bindings = expr[1]
        body = expr[2]

        # Evaluate all bindings in the *current* scope first, independently,
        # so they don't see each other's local assignments (parallel binding).
        evaluated_bindings = {}
        for binding in bindings:
            var_name = binding[0]
            evaluated_bindings[var_name] = self.transform(binding[1])

        shadowed = {}
        for var_name, parsed_val in evaluated_bindings.items():
            if var_name in self._symbols:
                shadowed[var_name] = self._symbols[var_name]
            self._symbols[var_name] = parsed_val

        result = self.transform(body)

        for var_name in evaluated_bindings.keys():
            if var_name in shadowed:
                self._symbols[var_name] = shadowed[var_name]
            else:
                del self._symbols[var_name]

        return result

    def transform_forall(self, expr):
        """
        Generate a wrapped SymPy ForAll AST node over bounds.

        Parameters
        ==========

        expr : list
            A list containing the command, scoped variables, and body

        Returns
        =======

        result : Function
            An uninterpreted function node evaluating the quantifier
        """
        bindings = expr[1]
        body = expr[2]

        shadowed = {}
        var_symbols = []
        for binding in bindings:
            var_name = binding[0]
            var_sym = Symbol(var_name)
            var_symbols.append(var_sym)
            if var_name in self._symbols:
                shadowed[var_name] = self._symbols[var_name]
            self._symbols[var_name] = var_sym

        parsed_body = self.transform(body)

        for binding in bindings:
            var_name = binding[0]
            if var_name in shadowed:
                self._symbols[var_name] = shadowed[var_name]
            else:
                del self._symbols[var_name]

        # Wrap the parsed body in a symbolic 'forall' function node, as
        # SymPy lacks native first-order logic quantifier nodes.
        return Function('forall')(tuple(var_symbols), parsed_body)

    def transform_exists(self, expr):
        """
        Generate a wrapped SymPy Exists AST node over bounds.

        Parameters
        ==========

        expr : list
            A list containing the command, scoped variables, and body

        Returns
        =======

        result : Function
            An uninterpreted function node evaluating the quantifier
        """
        bindings = expr[1]
        body = expr[2]

        shadowed = {}
        var_symbols = []
        for binding in bindings:
            var_name = binding[0]
            var_sym = Symbol(var_name)
            var_symbols.append(var_sym)
            if var_name in self._symbols:
                shadowed[var_name] = self._symbols[var_name]
            self._symbols[var_name] = var_sym

        parsed_body = self.transform(body)

        for binding in bindings:
            var_name = binding[0]
            if var_name in shadowed:
                self._symbols[var_name] = shadowed[var_name]
            else:
                del self._symbols[var_name]

        return Function('exists')(tuple(var_symbols), parsed_body)

    def transform_as(self, expr):
        """
        Bypass the type casting annotation natively in SymPy logic.

        Parameters
        ==========

        expr : list
            A list containing the command, expression, and sort

        Returns
        =======

        result : Any
            The recursively parsed expression without the cast
        """
        # SymPy handles typing via Symbol assumptions (e.g., real=True).
        # We can safely bypass explicit casting and evaluate the raw term.
        return self.transform(expr[1])

    def transform_operator(self, expr):
        """
        Route unhandled generic operators to specific SymPy classes.

        Parameters
        ==========

        expr : list
            A list containing the operator and its arguments

        Returns
        =======

        result : Any
            The fully populated SymPy operational AST
        """
        op = expr[0]

        if op == '!':
            return self.transform(expr[1])

        args = [self.transform(arg) for arg in expr[1:]]

        if op == '=':
            if len(args) == 2:
                return Eq(args[0], args[1])
            return And(
                *(Eq(args[i], args[i+1]) for i in range(len(args)-1))
            )
        if op == 'not':
            return Not(*args)
        if op == 'and':
            return And(*args)
        if op == 'or':
            return Or(*args)
        if op == '=>':
            # Implication is strictly right-associative when chained.
            res = args[-1]
            for arg in reversed(args[:-1]):
                res = Implies(arg, res)
            return res
        if op == '<->':
            if len(args) == 2:
                return Equivalent(args[0], args[1])
            # For >2 args, <-> means pairwise equivalence between adjacent
            # elements
            return And(
                *(Equivalent(args[i], args[i+1]) for i in range(len(args)-1))
            )
        if op == 'xor':
            return Xor(*args)
        if op == 'distinct':
            if len(args) == 2:
                return Ne(args[0], args[1])
            # Multi-arg distinct means every element is mutually unequal.
            # We expand this into an And over all pairwise inequalities.
            return And(
                *(Ne(args[i], args[j]) for i in range(len(args))
                  for j in range(i+1, len(args)))
            )

        if op == 'ite':
            condition, true_val, false_val = args
            return Piecewise((true_val, condition), (false_val, True))

        if op == '+':
            return Add(*args)
        if op == '-':
            if len(args) == 1:
                # Unary negation
                return -args[0]
            # Left-associative subtraction
            return args[0] - Add(*args[1:])
        if op == '*':
            return Mul(*args)
        if op == '/':
            if len(args) == 2:
                return args[0] / args[1]
            raise NotImplementedError
        if op == 'abs':
            return Abs(args[0])

        if op == '<':
            if len(args) == 2:
                return StrictLessThan(args[0], args[1])
            return And(
                *(StrictLessThan(args[i], args[i+1])
                  for i in range(len(args)-1))
            )
        if op == '>':
            if len(args) == 2:
                return StrictGreaterThan(args[0], args[1])
            return And(
                *(StrictGreaterThan(args[i], args[i+1])
                  for i in range(len(args)-1))
            )
        if op == '<=':
            if len(args) == 2:
                return LessThan(args[0], args[1])
            return And(
                *(LessThan(args[i], args[i+1]) for i in range(len(args)-1))
            )
        if op == '>=':
            if len(args) == 2:
                return GreaterThan(args[0], args[1])
            return And(
                *(GreaterThan(args[i], args[i+1]) for i in range(len(args)-1))
            )

        if op in self._macros:
            arg_names, temp_syms, parsed_body = self._macros[op]
            sub_dict = {
                temp_syms[name]: val for name, val in zip(arg_names, args)
            }
            return parsed_body.subs(sub_dict)

        if op in self._functions:
            return self._functions[op](*args)

        if op in self._extended_ops:
            # Wrap unknown but valid extended theory operations (e.g.
            # BitVectors)
            # in uninterpreted SymPy Function nodes to preserve the AST
            # structure.
            return Function(op)(*args)

        raise UnknownSMTLibCommandError(f"Operator '{op}' not supported.")


def parse_smtlib(source):
    """
    Main entrypoint wrapper to parse an SMT-LIB formatted string.

    Parameters
    ==========

    source : str
        A string containing the SMT-LIB source code

    Returns
    =======

    (symbols, assertions)
        A tuple of the parsed global symbol map and assertion tree
    """
    parser = SMTLibParser()
    return parser.parse_str(source)
