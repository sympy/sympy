"""Options manager for :class:`Poly` and public API functions. """

from sympy.core import S, Basic, sympify
from sympy.utilities import numbered_symbols, topological_sort
from sympy.utilities.iterables import has_dups

from sympy.polys.polyerrors import (
    GeneratorsError,
    OptionError,
    FlagError,
)

import sympy.polys

import re

class Option(object):
    """Base class for all kinds of options. """

    option = None

    is_Flag = False

    requires = []
    excludes = []

    after = []
    before = []

    @classmethod
    def default(cls):
        return None

    @classmethod
    def preprocess(cls, option):
        return None

    @classmethod
    def postprocess(cls, options):
        pass

class Flag(Option):
    """Base class for all kinds of flags. """

    is_Flag = True

class BooleanOption(Option):
    """An option that must have a boolean value or equivalent assigned. """

    @classmethod
    def preprocess(cls, value):
        if value in [True, False]:
            return bool(value)
        else:
            raise OptionError("'%s' must have a boolean value assigned, got %s" % (cls.option, value))

class OptionType(type):
    """Base type for all options that does registers options. """

    def __init__(cls, *args, **kwargs):
        @property
        def getter(self):
            try:
                return self[cls.option]
            except KeyError:
                return cls.default()

        setattr(Options, cls.option, getter)
        Options.__options__[cls.option] = cls

class Options(dict):
    """
    Options manager for polynomial manipulation module.

    Examples
    ========

    >>> from sympy.polys.polyoptions import Options
    >>> from sympy.polys.polyoptions import build_options

    >>> from sympy.abc import x, y, z

    >>> Options((x, y, z), {'domain': 'ZZ'})
    {'auto': False, 'domain': ZZ, 'gens': (x, y, z)}

    >>> build_options((x, y, z), {'domain': 'ZZ'})
    {'auto': False, 'domain': ZZ, 'gens': (x, y, z)}

    **Options**

    * Expand --- boolean option
    * Gens --- option
    * Wrt --- option
    * Sort --- option
    * Order --- option
    * Field --- boolean option
    * Greedy --- boolean option
    * Domain --- option
    * Split --- boolean option
    * Gaussian --- boolean option
    * Extension --- option
    * Modulus --- option
    * Symmetric --- boolean option
    * Strict --- boolean option

    **Flags**

    * Auto --- boolean flag
    * Frac --- boolean flag
    * Formal --- boolean flag
    * Polys --- boolean flag
    * Include --- boolean flag
    * All --- boolean flag
    * Gen --- flag

    """

    __order__ = None
    __options__ = {}

    def __init__(self, gens, args, flags=None, strict=False):
        dict.__init__(self)

        if gens and args.get('gens', ()):
            raise OptionError("both '*gens' and keyword argument 'gens' supplied")
        elif gens:
            args = dict(args)
            args['gens'] = gens

        defaults = args.pop('defaults', {})

        def preprocess_options(args):
            for option, value in args.iteritems():
                try:
                    cls = self.__options__[option]
                except KeyError:
                    raise OptionError("'%s' is not a valid option" % option)

                if issubclass(cls, Flag):
                    if flags is None or option not in flags:
                        if strict:
                            raise OptionError("'%s' flag is not allowed in this context" % option)

                if value is not None:
                    self[option] = cls.preprocess(value)

        preprocess_options(args)

        for key, value in dict(defaults).iteritems():
            if key in self:
                del defaults[key]
            else:
                for option in self.keys():
                    cls = self.__options__[option]

                    if key in cls.excludes:
                        del defaults[key]
                        break

        preprocess_options(defaults)

        for option in self.keys():
            cls = self.__options__[option]

            for require_option in cls.requires:
                if self.get(require_option) is None:
                    raise OptionError("'%s' option is only allowed together with '%s'" % (option, require_option))

            for exclude_option in cls.excludes:
                if self.get(exclude_option) is not None:
                    raise OptionError("'%s' option is not allowed together with '%s'" % (option, exclude_option))

        for option in self.__order__:
            self.__options__[option].postprocess(self)

    @classmethod
    def _init_dependencies_order(cls):
        """Resolve the order of options' processing. """
        if cls.__order__ is None:
            vertices, edges = [], set([])

            for name, option in cls.__options__.iteritems():
                vertices.append(name)

                for _name in option.after:
                    edges.add((_name, name))

                for _name in option.before:
                    edges.add((name, _name))

            try:
                cls.__order__ = topological_sort((vertices, list(edges)))
            except ValueError:
                raise RuntimeError("cycle detected in sympy.polys options framework")

    def clone(self, updates={}):
        """Clone ``self`` and update specified options. """
        obj = dict.__new__(self.__class__)

        for option, value in self.iteritems():
            obj[option] = value

        for option, value in updates.iteritems():
            obj[option] = value

        return obj

    def __setattr__(self, attr, value):
        if attr in self.__options__:
            self[attr] = value
        else:
            super(Options, self).__setattr__(attr, value)

    @property
    def args(self):
        args = {}

        for option, value in self.iteritems():
            if value is not None and option != 'gens':
                cls = self.__options__[option]

                if not issubclass(cls, Flag):
                    args[option] = value

        return args

    @property
    def options(self):
        options = {}

        for option, cls in self.__options__.iteritems():
            if not issubclass(cls, Flag):
                options[option] = getattr(self, option)

        return options

    @property
    def flags(self):
        flags = {}

        for option, cls in self.__options__.iteritems():
            if issubclass(cls, Flag):
                flags[option] = getattr(self, option)

        return flags

class Expand(BooleanOption):
    """``expand`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'expand'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return True

class Gens(Option):
    """``gens`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'gens'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return ()

    @classmethod
    def preprocess(cls, gens):
        if isinstance(gens, Basic):
            gens = (gens,)
        elif len(gens) == 1 and hasattr(gens[0], '__iter__'):
            gens = gens[0]

        if gens == (None,):
            gens = ()
        elif has_dups(gens):
            raise GeneratorsError("duplicated generators: %s" % str(gens))
        elif any(gen.is_commutative is False for gen in gens):
            raise GeneratorsError("non-commutative generators: %s" % str(gens))

        return tuple(gens)

class Wrt(Option):
    """``wrt`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'wrt'

    requires = []
    excludes = []

    _re_split = re.compile(r"\s*,\s*|\s+")

    @classmethod
    def preprocess(cls, wrt):
        if isinstance(wrt, Basic):
            return [str(wrt)]
        elif isinstance(wrt, str):
            wrt = wrt.strip()
            if wrt.endswith(','):
                raise OptionError('Bad input: missing parameter.')
            if not wrt:
                return []
            return [ gen for gen in cls._re_split.split(wrt) ]
        elif hasattr(wrt, '__getitem__'):
            return list(map(str, wrt))
        else:
            raise OptionError("invalid argument for 'wrt' option")

class Sort(Option):
    """``sort`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'sort'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return []

    @classmethod
    def preprocess(cls, sort):
        if isinstance(sort, str):
            return [ gen.strip() for gen in sort.split('>') ]
        elif hasattr(sort, '__getitem__'):
            return list(map(str, sort))
        else:
            raise OptionError("invalid argument for 'sort' option")

class Order(Option):
    """``order`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'order'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return sympy.polys.monomialtools.lex

    @classmethod
    def preprocess(cls, order):
        return sympy.polys.monomialtools.monomial_key(order)

class Field(BooleanOption):
    """``field`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'field'

    requires = []
    excludes = ['domain', 'split', 'gaussian']

class Greedy(BooleanOption):
    """``greedy`` option to polynomial manipulation functions. """
    __metaclass__ = OptionType

    option = 'greedy'

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus', 'symmetric']

class Composite(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'composite'

    @classmethod
    def default(cls):
        return True

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus', 'symmetric']

class Domain(Option):
    """``domain`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'domain'

    requires = []
    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    _re_finitefield = re.compile("^(FF|GF)\((\d+)\)$")
    _re_polynomial  = re.compile("^(Z|ZZ|Q|QQ)\[(.+)\]$")
    _re_fraction    = re.compile("^(Z|ZZ|Q|QQ)\((.+)\)$")
    _re_algebraic   = re.compile("^(Q|QQ)\<(.+)\>$")

    @classmethod
    def preprocess(cls, domain):
        if not isinstance(domain, str):
            return domain
        else:
            if domain in ['Z', 'ZZ']:
                return sympy.polys.domains.ZZ

            if domain in ['Q', 'QQ']:
                return sympy.polys.domains.QQ

            if domain in ['R', 'RR']:
                return sympy.polys.domains.RR

            if domain == 'EX':
                return sympy.polys.domains.EX

            r = cls._re_finitefield.match(domain)

            if r is not None:
                return sympy.polys.domains.FF(int(r.groups()[1]))

            r = cls._re_polynomial.match(domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return sympy.polys.domains.ZZ.poly_ring(*gens)
                else:
                    return sympy.polys.domains.QQ.poly_ring(*gens)

            r = cls._re_fraction.match(domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return sympy.polys.domains.ZZ.frac_field(*gens)
                else:
                    return sympy.polys.domains.QQ.frac_field(*gens)

            r = cls._re_algebraic.match(domain)

            if r is not None:
                gens = map(sympify, r.groups()[1].split(','))
                return sympy.polys.domains.QQ.algebraic_field(*gens)

            raise OptionError('expected a valid domain specification, got %s' % domain)

    @classmethod
    def postprocess(cls, options):
        if 'gens' in options and 'domain' in options and options['domain'].is_Composite and \
                (set(options['domain'].gens) & set(options['gens'])):
            raise GeneratorsError("ground domain and generators interferes together")

class Split(BooleanOption):
    """``split`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'split'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'gaussian', 'extension', 'modulus', 'symmetric']

    @classmethod
    def postprocess(cls, options):
        if 'split' in options:
            raise NotImplementedError("'split' option is not implemented yet")

class Gaussian(BooleanOption):
    """``gaussian`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'gaussian'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'split', 'extension', 'modulus', 'symmetric']

    @classmethod
    def postprocess(cls, options):
        if 'gaussian' in options and options['gaussian'] is True:
            options['extension'] = set([S.ImaginaryUnit])
            Extension.postprocess(options)

class Extension(Option):
    """``extension`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'extension'

    requires = []
    excludes = ['greedy', 'domain', 'split', 'gaussian', 'modulus', 'symmetric']

    @classmethod
    def preprocess(cls, extension):
        if extension == 1:
            return bool(extension)
        elif extension == 0:
            raise OptionError("'False' is an invalid argument for 'extension'")
        else:
            if not hasattr(extension, '__iter__'):
                extension = set([extension])
            else:
                if not extension:
                    extension = None
                else:
                    extension = set(extension)

            return extension

    @classmethod
    def postprocess(cls, options):
        if 'extension' in options and options['extension'] is not True:
            options['domain'] = sympy.polys.domains.QQ.algebraic_field(*options['extension'])

class Modulus(Option):
    """``modulus`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'modulus'

    requires = []
    excludes = ['greedy', 'split', 'domain', 'gaussian', 'extension']

    @classmethod
    def preprocess(cls, modulus):
        modulus = sympify(modulus)

        if modulus.is_Integer and modulus > 0:
            return int(modulus)
        else:
            raise OptionError("'modulus' must a positive integer, got %s" % modulus)

    @classmethod
    def postprocess(cls, options):
        if 'modulus' in options:
            modulus = options['modulus']
            symmetric = options.get('symmetric', True)
            options['domain'] = sympy.polys.domains.FF(modulus, symmetric)

class Symmetric(BooleanOption):
    """``symmetric`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'symmetric'

    requires = ['modulus']
    excludes = ['greedy', 'domain', 'split', 'gaussian', 'extension']

class Strict(BooleanOption):
    """``strict`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'strict'

    @classmethod
    def default(cls):
        return True

class Auto(BooleanOption, Flag):
    """``auto`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'auto'

    after = ['field', 'domain', 'extension', 'gaussian']

    @classmethod
    def default(cls):
        return True

    @classmethod
    def postprocess(cls, options):
        if ('domain' in options or 'field' in options) and 'auto' not in options:
            options['auto'] = False

class Frac(BooleanOption, Flag):
    """``auto`` option to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'frac'

    @classmethod
    def default(cls):
        return False

class Formal(BooleanOption, Flag):
    """``formal`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'formal'

    @classmethod
    def default(cls):
        return False

class Polys(BooleanOption, Flag):
    """``polys`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'polys'

class Include(BooleanOption, Flag):
    """``include`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'include'

    @classmethod
    def default(cls):
        return False

class All(BooleanOption, Flag):
    """``all`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'all'

    @classmethod
    def default(cls):
        return False

class Gen(Flag):
    """``gen`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'gen'

    @classmethod
    def default(cls):
        return 0

    @classmethod
    def preprocess(cls, gen):
        if isinstance(gen, (Basic, int)):
            return gen
        else:
            raise OptionError("invalid argument for 'gen' option")

class Symbols(Flag):
    """``symbols`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'symbols'

    @classmethod
    def default(cls):
        return numbered_symbols('s', start=1)

    @classmethod
    def preprocess(cls, symbols):
        if hasattr(symbols, '__iter__'):
            return iter(symbols)
        else:
            raise OptionError("expected an iterator or iterable container, got %s" % symbols)

class Method(Flag):
    """``method`` flag to polynomial manipulation functions. """

    __metaclass__ = OptionType

    option = 'method'

    @classmethod
    def preprocess(cls, method):
        if isinstance(method, str):
            return method.lower()
        else:
            raise OptionError("expected a string, got %s" % method)

def build_options(gens, args=None):
    """Construct options from keyword arguments or ... options. """
    if args is None:
        gens, args = (), gens

    if len(args) != 1 or 'opt' not in args or gens:
        return Options(gens, args)
    else:
        return args['opt']

def allowed_flags(args, flags):
    """
    Allow specified flags to be used in the given context.

    Examples
    ========

    >>> from sympy.polys.polyoptions import allowed_flags
    >>> from sympy.polys.domains import ZZ

    >>> allowed_flags({'domain': ZZ}, [])

    >>> allowed_flags({'domain': ZZ, 'frac': True}, [])
    Traceback (most recent call last):
    ...
    FlagError: 'frac' flag is not allowed in this context

    >>> allowed_flags({'domain': ZZ, 'frac': True}, ['frac'])

    """
    flags = set(flags)

    for arg in args.iterkeys():
        try:
            if Options.__options__[arg].is_Flag and not arg in flags:
                raise FlagError("'%s' flag is not allowed in this context" % arg)
        except KeyError:
            raise OptionError("'%s' is not a valid option" % arg)

def set_defaults(options, **defaults):
    """Update options with default values. """
    if 'defaults' not in options:
        options = dict(options)
        options['defaults'] = defaults

    return options

Options._init_dependencies_order()
